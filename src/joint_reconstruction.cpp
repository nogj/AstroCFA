#include "astrocfa/joint_reconstruction.hpp"

#include "astrocfa/demosaic.hpp"

#include <algorithm>
#include <array>
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

struct AxisStencilEntry {
  std::size_t index;
  double weight;
};

struct AxisStencil {
  std::array<AxisStencilEntry, 128> entries;
  std::size_t size = 0;

  [[nodiscard]] const AxisStencilEntry *begin() const { return entries.data(); }
  [[nodiscard]] const AxisStencilEntry *end() const {
    return entries.data() + size;
  }
  [[nodiscard]] bool empty() const { return size == 0; }
};

void add_axis_weight(AxisStencil &stencil, std::size_t index,
                     double weight) {
  for(std::size_t i = 0; i < stencil.size; ++i) {
    AxisStencilEntry &entry = stencil.entries[i];
    if(entry.index == index) {
      entry.weight += weight;
      return;
    }
  }
  if(stencil.size >= stencil.entries.size()) {
    throw std::invalid_argument("PSF stencil exceeds supported radius");
  }
  stencil.entries[stencil.size++] =
      AxisStencilEntry{.index = index, .weight = weight};
}

AxisStencil make_axis_stencil(double position, std::size_t extent,
                              double psf_sigma) {
  AxisStencil stencil;
  const int radius = psf_sigma > 0.0
                         ? std::max(1, static_cast<int>(std::ceil(3.0 * psf_sigma)))
                         : 0;
  double weight_sum = 0.0;
  for(int tap = -radius; tap <= radius; ++tap) {
    const double shifted = position + static_cast<double>(tap);
    if(shifted < 0.0 || shifted > static_cast<double>(extent - 1U)) {
      continue;
    }
    const double psf_weight = psf_sigma > 0.0
                                  ? std::exp(-0.5 * tap * tap /
                                             (psf_sigma * psf_sigma))
                                  : 1.0;
    const std::size_t lower = static_cast<std::size_t>(std::floor(shifted));
    const std::size_t upper = std::min(lower + 1U, extent - 1U);
    const double fraction = shifted - static_cast<double>(lower);
    add_axis_weight(stencil, lower, psf_weight * (1.0 - fraction));
    if(fraction > 0.0) {
      add_axis_weight(stencil, upper, psf_weight * fraction);
    }
    weight_sum += psf_weight;
  }
  if(weight_sum > 0.0) {
    for(std::size_t i = 0; i < stencil.size; ++i) {
      stencil.entries[i].weight /= weight_sum;
    }
  }
  return stencil;
}

template <typename Visitor>
bool visit_measurement_stencil(double x, double y, double psf_sigma,
                               std::size_t width, std::size_t height,
                               Visitor visitor) {
  if(psf_sigma <= 0.0) {
    return visit_bilinear(x, y, width, height, visitor);
  }
  const AxisStencil x_stencil = make_axis_stencil(x, width, psf_sigma);
  const AxisStencil y_stencil = make_axis_stencil(y, height, psf_sigma);
  if(x_stencil.empty() || y_stencil.empty()) {
    return false;
  }
  for(const AxisStencilEntry &y_entry : y_stencil) {
    for(const AxisStencilEntry &x_entry : x_stencil) {
      visitor(x_entry.index, y_entry.index, x_entry.weight * y_entry.weight);
    }
  }
  return true;
}

double sample_measurement(const astrocfa::RgbImage &image, int channel, double x,
                          double y, double psf_sigma) {
  double value = 0.0;
  const bool inside = visit_measurement_stencil(
      x, y, psf_sigma, image.width(), image.height(),
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
  double reduced_chi_square = 0.0;
  std::size_t samples = 0;
  std::size_t outliers = 0;
};

ResidualSummary measure_residual(const astrocfa::RgbImage &image,
                                 const std::vector<astrocfa::JointCfaFrame> &frames,
                                 const astrocfa::JointReconstructionOptions &options) {
  ResidualSummary summary;
  double squared = 0.0;
  double normalized = 0.0;
  double robust_chi_square = 0.0;
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
        const double residual =
            static_cast<double>(sample.value) -
            sample_measurement(image, channel, target_x, target_y,
                               input.psf_sigma * scale);
        const astrocfa::NoiseEstimate noise =
            astrocfa::estimate_noise(sample.value, options.noise);
        squared += residual * residual;
        const double normalized_residual = residual / noise.sigma;
        normalized += std::abs(normalized_residual);
        robust_chi_square +=
            std::min(normalized_residual * normalized_residual,
                     options.huber_sigma * options.huber_sigma);
        summary.outliers +=
            std::abs(residual) > options.huber_sigma * noise.sigma ? 1U : 0U;
        summary.samples += 1;
      }
    }
  }
  if(summary.samples > 0) {
    summary.rmse = std::sqrt(squared / static_cast<double>(summary.samples));
    summary.normalized_mae = normalized / static_cast<double>(summary.samples);
    summary.reduced_chi_square =
        robust_chi_square / static_cast<double>(summary.samples);
  }
  return summary;
}

double luminance(astrocfa::RgbPixel pixel) {
  return 0.2126 * pixel.r + 0.7152 * pixel.g + 0.0722 * pixel.b;
}

void regularize_flux_conserving_luma(
    astrocfa::RgbImage &image,
    const astrocfa::JointReconstructionOptions &options) {
  const std::size_t width = image.width();
  const std::size_t height = image.height();
  std::vector<astrocfa::RgbPixel> delta(width * height);
  const double luma_step = 0.25 * options.luma_smoothness;
  const auto diffuse_pair = [&](std::size_t first_x, std::size_t first_y,
                                std::size_t second_x, std::size_t second_y) {
    const astrocfa::RgbPixel first = image.pixel(first_x, first_y);
    const astrocfa::RgbPixel second = image.pixel(second_x, second_y);
    const double conductance = std::exp(
        -options.edge_sensitivity * std::abs(luminance(second) - luminance(first)));
    const double green_transfer =
        luma_step * conductance * (static_cast<double>(second.g) - first.g);
    astrocfa::RgbPixel &first_delta = delta[first_y * width + first_x];
    astrocfa::RgbPixel &second_delta = delta[second_y * width + second_x];
    first_delta.g += static_cast<float>(green_transfer);
    second_delta.g -= static_cast<float>(green_transfer);
  };

  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      if(x + 1U < width) {
        diffuse_pair(x, y, x + 1U, y);
      }
      if(y + 1U < height) {
        diffuse_pair(x, y, x, y + 1U);
      }
    }
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const astrocfa::RgbPixel current = image.pixel(x, y);
      const astrocfa::RgbPixel change = delta[y * width + x];
      const double green = static_cast<double>(current.g) + change.g;
      const double red_chroma = static_cast<double>(current.r) - current.g;
      const double blue_chroma = static_cast<double>(current.b) - current.g;
      image.set_pixel(
          x, y,
          astrocfa::RgbPixel{
              .r = static_cast<float>(std::clamp(green + red_chroma, 0.0, 1.25)),
              .g = static_cast<float>(std::clamp(green, 0.0, 1.25)),
              .b = static_cast<float>(std::clamp(green + blue_chroma, 0.0, 1.25)),
          });
    }
  }
}

void regularize_chroma(astrocfa::RgbImage &image,
                       const astrocfa::JointReconstructionOptions &options,
                       const std::vector<astrocfa::RgbPixel> &data_support,
                       bool preserve_direct_measurements, bool psf_aware) {
  if((options.chroma_smoothness <= 0.0 && options.luma_smoothness <= 0.0) ||
     image.width() < 3 || image.height() < 3) {
    return;
  }
  if(psf_aware) {
    regularize_flux_conserving_luma(image, options);
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
        const double luma_amount =
            psf_aware ? 0.0 : std::clamp(options.luma_smoothness, 0.0, 1.0);
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
  if(options.scale == 0 || options.scale > 2) {
    throw std::invalid_argument("Joint CFA scale must be 1 or 2");
  }
  if(options.learning_rate <= 0.0 || options.learning_rate > 1.0 ||
     options.huber_sigma <= 0.0 || options.luma_smoothness < 0.0 ||
     options.luma_smoothness > 1.0 || options.chroma_smoothness < 0.0 ||
     options.chroma_smoothness > 1.0 || !std::isfinite(options.huber_sigma) ||
     !std::isfinite(options.learning_rate) || options.edge_sensitivity < 0.0 ||
     !std::isfinite(options.edge_sensitivity) ||
     options.discrepancy_target <= 0.0 || options.discrepancy_tolerance < 0.0 ||
     !std::isfinite(options.discrepancy_target) ||
     !std::isfinite(options.discrepancy_tolerance) ||
     !std::isfinite(options.luma_smoothness) ||
     !std::isfinite(options.chroma_smoothness)) {
    throw std::invalid_argument("Invalid joint reconstruction solver options");
  }

  const CfaFrame &reference = *frames.front().cfa;
  for(const auto &input : frames) {
    if(input.cfa == nullptr || input.cfa->width() != reference.width() ||
       input.cfa->height() != reference.height() || input.weight <= 0.0 ||
       input.psf_sigma < 0.0 || input.psf_sigma > 8.0 ||
       !std::isfinite(input.psf_sigma) || !std::isfinite(input.weight) ||
       !std::isfinite(input.offset.dx) || !std::isfinite(input.offset.dy)) {
      throw std::invalid_argument(
          "Joint CFA inputs require equal dimensions, positive weights, and PSF sigma in 0..8");
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
  const bool psf_aware = std::any_of(
      frames.begin(), frames.end(),
      [](const JointCfaFrame &input) { return input.psf_sigma > 0.0; });

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
            const double prediction = sample_measurement(
                image, channel, target_x, target_y, input.psf_sigma * scale);
            value -= prediction;
            const NoiseEstimate noise = estimate_noise(sample.value, options.noise);
            const double normalized = std::abs(value) / noise.sigma;
            const double robust = normalized > options.huber_sigma
                                      ? options.huber_sigma / normalized
                                      : 1.0;
            measurement_weight *= robust * noise.weight;
          }
          const double operator_sigma =
              residual_mode ? input.psf_sigma * scale : 0.0;
          visit_measurement_stencil(
              target_x, target_y, operator_sigma, width, height,
              [&](std::size_t px, std::size_t py, double operator_weight) {
                const std::size_t index = py * width + px;
                const double weighted = measurement_weight * operator_weight;
                channel_ref(numerator[index], channel) +=
                    static_cast<float>(weighted * value);
                channel_ref(denominator[index], channel) +=
                    static_cast<float>(residual_mode && !psf_aware
                                           ? weighted * operator_weight
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
  std::size_t executed_iterations = 0;
  bool stopped_by_discrepancy = false;
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
    const bool preserve_direct_measurements =
        frames.size() == 1U && frames.front().psf_sigma == 0.0;
    regularize_chroma(image, options, denominator, preserve_direct_measurements,
                      psf_aware);
    executed_iterations = iteration + 1U;
    if(options.stop_on_discrepancy &&
       executed_iterations >= options.minimum_iterations) {
      const ResidualSummary current = measure_residual(image, frames, options);
      const double threshold = options.discrepancy_target *
                               (1.0 + options.discrepancy_tolerance);
      if(current.reduced_chi_square <= threshold) {
        stopped_by_discrepancy = true;
        break;
      }
    }
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
  stats.iterations = executed_iterations;
  stats.maximum_iterations = options.iterations;
  stats.robust_outliers = final.outliers;
  for(const JointCfaFrame &input : frames) {
    if(input.psf_sigma <= 0.0) {
      continue;
    }
    if(stats.psf_frames == 0) {
      stats.minimum_psf_sigma = input.psf_sigma;
    } else {
      stats.minimum_psf_sigma =
          std::min(stats.minimum_psf_sigma, input.psf_sigma);
    }
    stats.maximum_psf_sigma =
        std::max(stats.maximum_psf_sigma, input.psf_sigma);
    stats.psf_frames += 1;
  }
  stats.initial_rmse = initial.rmse;
  stats.final_rmse = final.rmse;
  stats.final_normalized_mae = final.normalized_mae;
  stats.initial_reduced_chi_square = initial.reduced_chi_square;
  stats.final_reduced_chi_square = final.reduced_chi_square;
  stats.stopped_by_discrepancy = stopped_by_discrepancy;
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
