#include "astrocfa/joint_reconstruction.hpp"

#include "astrocfa/demosaic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

int channel_index(astrocfa::CfaColor color) {
  if(color == astrocfa::CfaColor::red) {
    return 0;
  }
  if(color == astrocfa::CfaColor::blue) {
    return 2;
  }
  return 1;
}

float get_channel(astrocfa::RgbPixel pixel, int channel) {
  return channel == 0 ? pixel.r : channel == 1 ? pixel.g : pixel.b;
}

float &channel_ref(astrocfa::RgbPixel &pixel, int channel) {
  return channel == 0 ? pixel.r : channel == 1 ? pixel.g : pixel.b;
}

bool compatible_pattern(const astrocfa::BayerPattern &first,
                        const astrocfa::BayerPattern &second) {
  for(std::size_t y = 0; y < 2; ++y) {
    for(std::size_t x = 0; x < 2; ++x) {
      if(first.at(x, y) != second.at(x, y)) {
        return false;
      }
    }
  }
  return true;
}

template <typename Visitor>
bool visit_bilinear(double x, double y, std::size_t width, std::size_t height,
                    Visitor visitor) {
  if(x < 0.0 || y < 0.0 || x > static_cast<double>(width - 1U) ||
     y > static_cast<double>(height - 1U)) {
    return false;
  }
  const std::size_t x0 = static_cast<std::size_t>(std::floor(x));
  const std::size_t y0 = static_cast<std::size_t>(std::floor(y));
  const std::size_t x1 = std::min(x0 + 1U, width - 1U);
  const std::size_t y1 = std::min(y0 + 1U, height - 1U);
  const double fx = x - static_cast<double>(x0);
  const double fy = y - static_cast<double>(y0);
  const double weights[4] = {
      (1.0 - fx) * (1.0 - fy), fx * (1.0 - fy),
      (1.0 - fx) * fy, fx * fy,
  };
  const std::size_t xs[4] = {x0, x1, x0, x1};
  const std::size_t ys[4] = {y0, y0, y1, y1};
  for(int i = 0; i < 4; ++i) {
    if(weights[i] > 0.0) {
      visitor(xs[i], ys[i], weights[i]);
    }
  }
  return true;
}

double sample_channel(const astrocfa::RgbImage &image, int channel, double x,
                      double y) {
  double value = 0.0;
  const bool inside = visit_bilinear(
      x, y, image.width(), image.height(),
      [&](std::size_t px, std::size_t py, double weight) {
        value += weight * get_channel(image.pixel(px, py), channel);
      });
  return inside ? value : 0.0;
}

astrocfa::RgbPixel bilinear_rgb(const astrocfa::RgbImage &image, double x,
                               double y) {
  astrocfa::RgbPixel value;
  visit_bilinear(x, y, image.width(), image.height(),
                 [&](std::size_t px, std::size_t py, double weight) {
                   const astrocfa::RgbPixel source = image.pixel(px, py);
                   value.r += static_cast<float>(weight * source.r);
                   value.g += static_cast<float>(weight * source.g);
                   value.b += static_cast<float>(weight * source.b);
                 });
  return value;
}

struct ResidualSummary {
  double rmse = 0.0;
  double normalized_mae = 0.0;
  std::size_t samples = 0;
  std::size_t outliers = 0;
};

ResidualSummary measure_residual(const astrocfa::RgbImage &image,
                                 const std::vector<astrocfa::JointCfaFrame> &frames,
                                 const astrocfa::JointReconstructionOptions &options) {
  ResidualSummary summary;
  double squared = 0.0;
  double normalized = 0.0;
  const double scale = static_cast<double>(options.scale);
  for(const auto &input : frames) {
    for(std::size_t y = 0; y < input.cfa->height(); ++y) {
      for(std::size_t x = 0; x < input.cfa->width(); ++x) {
        const astrocfa::CfaSample sample = input.cfa->sample_info(x, y);
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const double target_x = (static_cast<double>(x) + input.offset.dx) * scale;
        const double target_y = (static_cast<double>(y) + input.offset.dy) * scale;
        if(target_x < 0.0 || target_y < 0.0 ||
           target_x > static_cast<double>(image.width() - 1U) ||
           target_y > static_cast<double>(image.height() - 1U)) {
          continue;
        }
        const int channel = channel_index(input.cfa->pattern().at(x, y));
        const double residual = static_cast<double>(sample.value) -
                                sample_channel(image, channel, target_x, target_y);
        const astrocfa::NoiseEstimate noise =
            astrocfa::estimate_noise(sample.value, options.noise);
        squared += residual * residual;
        normalized += std::abs(residual) / noise.sigma;
        summary.outliers +=
            std::abs(residual) > options.huber_sigma * noise.sigma ? 1U : 0U;
        summary.samples += 1;
      }
    }
  }
  if(summary.samples > 0) {
    summary.rmse = std::sqrt(squared / static_cast<double>(summary.samples));
    summary.normalized_mae = normalized / static_cast<double>(summary.samples);
  }
  return summary;
}

double luminance(astrocfa::RgbPixel pixel) {
  return 0.2126 * pixel.r + 0.7152 * pixel.g + 0.0722 * pixel.b;
}

void regularize_chroma(astrocfa::RgbImage &image,
                       const astrocfa::JointReconstructionOptions &options,
                       const std::vector<astrocfa::RgbPixel> &data_support,
                       bool preserve_direct_measurements) {
  if((options.chroma_smoothness <= 0.0 && options.luma_smoothness <= 0.0) ||
     image.width() < 3 || image.height() < 3) {
    return;
  }
  const int dx[4] = {-1, 1, 0, 0};
  const int dy[4] = {0, 0, -1, 1};
  for(int parity = 0; parity < 2; ++parity) {
    for(std::size_t y = 1; y + 1 < image.height(); ++y) {
      for(std::size_t x = 1; x + 1 < image.width(); ++x) {
        if(static_cast<int>((x + y) & 1U) != parity) {
          continue;
        }
        astrocfa::RgbPixel center = image.pixel(x, y);
        const double center_luma = luminance(center);
        double red_chroma = 0.0;
        double blue_chroma = 0.0;
        double green = 0.0;
        double weight_sum = 0.0;
        for(int i = 0; i < 4; ++i) {
          const astrocfa::RgbPixel neighbor = image.pixel(
              static_cast<std::size_t>(static_cast<int>(x) + dx[i]),
              static_cast<std::size_t>(static_cast<int>(y) + dy[i]));
          const double weight = std::exp(
              -options.edge_sensitivity * std::abs(luminance(neighbor) - center_luma));
          red_chroma += weight * (neighbor.r - neighbor.g);
          blue_chroma += weight * (neighbor.b - neighbor.g);
          green += weight * neighbor.g;
          weight_sum += weight;
        }
        if(weight_sum <= 0.0) {
          continue;
        }
        const double chroma_amount = std::clamp(options.chroma_smoothness, 0.0, 1.0);
        const double luma_amount = std::clamp(options.luma_smoothness, 0.0, 1.0);
        const double current_red = center.r - center.g;
        const double current_blue = center.b - center.g;
        const astrocfa::RgbPixel support = data_support[y * image.width() + x];
        const double regularized_green =
            (!preserve_direct_measurements || support.g <= 0.0F)
                ? (1.0 - luma_amount) * center.g + luma_amount * green / weight_sum
                : center.g;
        if(!preserve_direct_measurements || support.r <= 0.0F) {
          center.r = static_cast<float>(std::clamp(
              regularized_green + (1.0 - chroma_amount) * current_red +
                  chroma_amount * red_chroma / weight_sum,
              0.0, 1.25));
        }
        if(!preserve_direct_measurements || support.b <= 0.0F) {
          center.b = static_cast<float>(std::clamp(
              regularized_green + (1.0 - chroma_amount) * current_blue +
                  chroma_amount * blue_chroma / weight_sum,
              0.0, 1.25));
        }
        center.g = static_cast<float>(std::clamp(regularized_green, 0.0, 1.25));
        image.set_pixel(x, y, center);
      }
    }
  }
}

} // namespace

namespace astrocfa {

JointReconstructionResult reconstruct_joint_cfa(
    const std::vector<JointCfaFrame> &frames, JointReconstructionOptions options) {
  if(frames.empty() || frames.front().cfa == nullptr) {
    throw std::invalid_argument("Joint CFA reconstruction requires input frames");
  }
  if(options.scale == 0) {
    throw std::invalid_argument("Joint CFA scale must be positive");
  }
  if(options.learning_rate <= 0.0 || options.learning_rate > 1.0 ||
     options.huber_sigma <= 0.0 || options.luma_smoothness < 0.0 ||
     options.luma_smoothness > 1.0 || options.chroma_smoothness < 0.0 ||
     options.chroma_smoothness > 1.0 || !std::isfinite(options.huber_sigma) ||
     !std::isfinite(options.luma_smoothness) ||
     !std::isfinite(options.chroma_smoothness)) {
    throw std::invalid_argument("Invalid joint reconstruction solver options");
  }

  const CfaFrame &reference = *frames.front().cfa;
  for(const auto &input : frames) {
    if(input.cfa == nullptr || input.cfa->width() != reference.width() ||
       input.cfa->height() != reference.height() || input.weight <= 0.0) {
      throw std::invalid_argument("Joint CFA inputs must be non-null, positive-weight, and equal-sized");
    }
    if(!compatible_pattern(reference.pattern(), input.cfa->pattern())) {
      throw std::invalid_argument("Joint CFA input Bayer phases differ");
    }
  }

  const std::size_t width = reference.width() * options.scale;
  const std::size_t height = reference.height() * options.scale;
  const std::size_t pixels = width * height;
  RgbImage image(width, height);
  {
    const RgbImage baseline = reconstruct_frequency_guided(reference, options.noise).image;
    for(std::size_t y = 0; y < height; ++y) {
      for(std::size_t x = 0; x < width; ++x) {
        image.set_pixel(x, y,
                        bilinear_rgb(baseline,
                                     static_cast<double>(x) / options.scale,
                                     static_cast<double>(y) / options.scale));
      }
    }
  }

  std::vector<RgbPixel> numerator(pixels);
  std::vector<RgbPixel> denominator(pixels);
  const auto clear_accumulators = [&]() {
    std::fill(numerator.begin(), numerator.end(), RgbPixel{});
    std::fill(denominator.begin(), denominator.end(), RgbPixel{});
  };
  const double scale = static_cast<double>(options.scale);

  const auto accumulate_measurements = [&](bool residual_mode) {
    clear_accumulators();
    for(const auto &input : frames) {
      for(std::size_t y = 0; y < input.cfa->height(); ++y) {
        for(std::size_t x = 0; x < input.cfa->width(); ++x) {
          const CfaSample sample = input.cfa->sample_info(x, y);
          if(!sample.valid || sample.clipped) {
            continue;
          }
          const double target_x = (static_cast<double>(x) + input.offset.dx) * scale;
          const double target_y = (static_cast<double>(y) + input.offset.dy) * scale;
          const int channel = channel_index(input.cfa->pattern().at(x, y));
          double value = sample.value;
          double measurement_weight = input.weight;
          if(residual_mode) {
            const double prediction = sample_channel(image, channel, target_x, target_y);
            value -= prediction;
            const NoiseEstimate noise = estimate_noise(sample.value, options.noise);
            const double normalized = std::abs(value) / noise.sigma;
            const double robust = normalized > options.huber_sigma
                                      ? options.huber_sigma / normalized
                                      : 1.0;
            measurement_weight *= robust * noise.weight;
          }
          visit_bilinear(target_x, target_y, width, height,
                         [&](std::size_t px, std::size_t py, double kernel_weight) {
                           const std::size_t index = py * width + px;
                           const double weighted = measurement_weight * kernel_weight;
                           channel_ref(numerator[index], channel) +=
                               static_cast<float>(weighted * value);
                           channel_ref(denominator[index], channel) +=
                               static_cast<float>(residual_mode
                                                      ? weighted * kernel_weight
                                                      : weighted);
                         });
        }
      }
    }
  };

  accumulate_measurements(false);
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const std::size_t index = y * width + x;
      RgbPixel pixel = image.pixel(x, y);
      for(int channel = 0; channel < 3; ++channel) {
        const float weight = get_channel(denominator[index], channel);
        if(weight > 0.0F) {
          channel_ref(pixel, channel) = get_channel(numerator[index], channel) / weight;
        }
      }
      image.set_pixel(x, y, pixel);
    }
  }

  const ResidualSummary initial = measure_residual(image, frames, options);
  for(std::size_t iteration = 0; iteration < options.iterations; ++iteration) {
    accumulate_measurements(true);
    for(std::size_t y = 0; y < height; ++y) {
      for(std::size_t x = 0; x < width; ++x) {
        const std::size_t index = y * width + x;
        RgbPixel pixel = image.pixel(x, y);
        for(int channel = 0; channel < 3; ++channel) {
          const double normalizer = get_channel(denominator[index], channel);
          if(normalizer > 0.0) {
            const double delta = get_channel(numerator[index], channel) / normalizer;
            channel_ref(pixel, channel) = static_cast<float>(std::clamp(
                static_cast<double>(get_channel(pixel, channel)) +
                    options.learning_rate * delta,
                0.0, 1.25));
          }
        }
        image.set_pixel(x, y, pixel);
      }
    }
    regularize_chroma(image, options, denominator, frames.size() == 1U);
  }

  accumulate_measurements(true);
  std::array<double, 3> max_weight = {0.0, 0.0, 0.0};
  std::array<std::size_t, 3> covered = {0, 0, 0};
  for(const RgbPixel &weight : denominator) {
    for(int channel = 0; channel < 3; ++channel) {
      const double value = get_channel(weight, channel);
      max_weight[channel] = std::max(max_weight[channel], value);
      covered[channel] += value > 0.0 ? 1U : 0U;
    }
  }
  numerator.clear();
  numerator.shrink_to_fit();
  RgbImage confidence(width, height);
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const RgbPixel weight = denominator[y * width + x];
      confidence.set_pixel(
          x, y,
          RgbPixel{
              .r = static_cast<float>(max_weight[0] > 0.0 ? weight.r / max_weight[0] : 0.0),
              .g = static_cast<float>(max_weight[1] > 0.0 ? weight.g / max_weight[1] : 0.0),
              .b = static_cast<float>(max_weight[2] > 0.0 ? weight.b / max_weight[2] : 0.0),
          });
    }
  }

  const ResidualSummary final = measure_residual(image, frames, options);
  JointReconstructionStats stats;
  stats.frames = frames.size();
  stats.measurements = final.samples;
  stats.iterations = options.iterations;
  stats.robust_outliers = final.outliers;
  stats.initial_rmse = initial.rmse;
  stats.final_rmse = final.rmse;
  stats.final_normalized_mae = final.normalized_mae;
  for(int channel = 0; channel < 3; ++channel) {
    stats.channel_coverage[channel] =
        pixels > 0 ? static_cast<double>(covered[channel]) / static_cast<double>(pixels)
                   : 0.0;
  }
  return JointReconstructionResult{
      .image = std::move(image),
      .confidence = std::move(confidence),
      .stats = stats,
  };
}

} // namespace astrocfa
