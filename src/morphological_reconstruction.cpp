#include "astrocfa/morphological_reconstruction.hpp"

#include "astrocfa/luminance_proxy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

struct SourceCandidate {
  double x = 0.0;
  double y = 0.0;
  double sigma = 0.7;
  double peak = 0.0;
  double evidence = 0.0;
  std::array<double, 3> amplitude{};
};

struct RadialProfile {
  double bin_size = 0.125;
  std::vector<double> values;
};

struct EpsfProfile {
  double spacing = 0.25;
  double radius = 4.0;
  std::size_t size = 0;
  std::vector<double> values;
};

struct ChromaticEpsfTransform {
  std::array<double, 3> shift_x{};
  std::array<double, 3> shift_y{};
  std::array<double, 3> scale{1.0, 1.0, 1.0};
};

struct ObservableSourceSelection {
  std::vector<SourceCandidate> sources;
  double information_gain = 0.0;
};

int channel_index(astrocfa::CfaColor color) {
  if(color == astrocfa::CfaColor::red) {
    return 0;
  }
  if(color == astrocfa::CfaColor::blue) {
    return 2;
  }
  return 1;
}

double channel_value(astrocfa::RgbPixel pixel, int channel) {
  return channel == 0 ? pixel.r : channel == 1 ? pixel.g : pixel.b;
}

void set_channel(astrocfa::RgbPixel &pixel, int channel, double value) {
  const float bounded = static_cast<float>(std::clamp(value, 0.0, 1.25));
  if(channel == 0) {
    pixel.r = bounded;
  } else if(channel == 1) {
    pixel.g = bounded;
  } else {
    pixel.b = bounded;
  }
}

std::vector<double> gaussian_kernel(double sigma) {
  const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));
  std::vector<double> kernel(static_cast<std::size_t>(2 * radius + 1));
  double sum = 0.0;
  for(int offset = -radius; offset <= radius; ++offset) {
    const double value = std::exp(-0.5 * offset * offset / (sigma * sigma));
    kernel[static_cast<std::size_t>(offset + radius)] = value;
    sum += value;
  }
  for(double &value : kernel) {
    value /= sum;
  }
  return kernel;
}

std::vector<double> gaussian_blur(const std::vector<double> &input,
                                  std::size_t width, std::size_t height,
                                  double sigma) {
  const std::vector<double> kernel = gaussian_kernel(sigma);
  const int radius = static_cast<int>(kernel.size() / 2U);
  std::vector<double> horizontal(input.size());
  std::vector<double> output(input.size());
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      double sum = 0.0;
      for(int offset = -radius; offset <= radius; ++offset) {
        const std::size_t sx = static_cast<std::size_t>(std::clamp(
            static_cast<int>(x) + offset, 0, static_cast<int>(width) - 1));
        sum += kernel[static_cast<std::size_t>(offset + radius)] *
               input[y * width + sx];
      }
      horizontal[y * width + x] = sum;
    }
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      double sum = 0.0;
      for(int offset = -radius; offset <= radius; ++offset) {
        const std::size_t sy = static_cast<std::size_t>(std::clamp(
            static_cast<int>(y) + offset, 0, static_cast<int>(height) - 1));
        sum += kernel[static_cast<std::size_t>(offset + radius)] *
               horizontal[sy * width + x];
      }
      output[y * width + x] = sum;
    }
  }
  return output;
}

bool validation_sample(int x, int y);

std::vector<SourceCandidate> detect_sources(
    const astrocfa::CfaFrame &cfa, const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options) {
  const std::size_t width = cfa.width() / 2U;
  const std::size_t height = cfa.height() / 2U;
  astrocfa::LuminanceProxy proxy(width, height);
  std::vector<float> samples;
  samples.reserve(width * height);
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      double sum = 0.0;
      std::size_t count = 0;
      for(std::size_t dy = 0; dy < 2; ++dy) {
        for(std::size_t dx = 0; dx < 2; ++dx) {
          const std::size_t sx = 2U * x + dx;
          const std::size_t sy = 2U * y + dy;
          if(validation_sample(static_cast<int>(sx), static_cast<int>(sy))) {
            continue;
          }
          const astrocfa::CfaSample sample = cfa.sample_info(sx, sy);
          if(sample.valid && !sample.clipped) {
            sum += sample.value;
            count += 1U;
          }
        }
      }
      const float value = count > 0 ? static_cast<float>(sum / count) : 0.0F;
      proxy.set_value(x, y, value);
      if(value > 0.0F) {
        samples.push_back(value);
      }
    }
  }
  const astrocfa::RobustBackground background =
      astrocfa::robust_background(samples);
  const double detector_sigma = astrocfa::estimate_noise(
                                    background.location, noise_model)
                                    .sigma * 0.5;
  const double threshold = background.location +
                           options.detection_sigma * detector_sigma;

  std::vector<SourceCandidate> candidates;
  for(std::size_t y = 1; y + 1 < height; ++y) {
    for(std::size_t x = 1; x + 1 < width; ++x) {
      const double value = proxy.value(x, y);
      if(value <= threshold) {
        continue;
      }
      bool local_maximum = true;
      for(int dy = -1; dy <= 1 && local_maximum; ++dy) {
        for(int dx = -1; dx <= 1; ++dx) {
          if(dx == 0 && dy == 0) {
            continue;
          }
          const std::size_t nx = static_cast<std::size_t>(static_cast<int>(x) + dx);
          const std::size_t ny = static_cast<std::size_t>(static_cast<int>(y) + dy);
          if(proxy.value(nx, ny) > value) {
            local_maximum = false;
            break;
          }
        }
      }
      if(local_maximum) {
        candidates.push_back(SourceCandidate{
            .x = 2.0 * static_cast<double>(x) + 0.5,
            .y = 2.0 * static_cast<double>(y) + 0.5,
            .peak = value - background.location,
        });
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const SourceCandidate &first, const SourceCandidate &second) {
              return first.peak > second.peak;
            });
  std::vector<SourceCandidate> selected;
  for(const SourceCandidate &candidate : candidates) {
    const bool overlaps = std::any_of(
        selected.begin(), selected.end(), [&](const SourceCandidate &source) {
          const double dx = source.x - candidate.x;
          const double dy = source.y - candidate.y;
          return dx * dx + dy * dy < 16.0;
        });
    if(!overlaps) {
      selected.push_back(candidate);
      if(selected.size() == options.maximum_sources) {
        break;
      }
    }
  }
  return selected;
}

bool inside_source_mask(std::size_t x, std::size_t y,
                        const std::vector<SourceCandidate> &sources,
                        double radius) {
  const double radius2 = radius * radius;
  return std::any_of(sources.begin(), sources.end(),
                     [&](const SourceCandidate &source) {
                       const double dx = static_cast<double>(x) - source.x;
                       const double dy = static_cast<double>(y) - source.y;
                       return dx * dx + dy * dy <= radius2;
                     });
}

astrocfa::RgbImage estimate_diffuse(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &fallback,
    const std::vector<SourceCandidate> &sources,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only = false) {
  const std::size_t width = cfa.width();
  const std::size_t height = cfa.height();
  std::array<std::vector<double>, 3> numerator;
  std::array<std::vector<double>, 3> denominator;
  for(int channel = 0; channel < 3; ++channel) {
    numerator[channel].assign(width * height, 0.0);
    denominator[channel].assign(width * height, 0.0);
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const astrocfa::CfaSample sample = cfa.sample_info(x, y);
      if((training_only &&
          validation_sample(static_cast<int>(x), static_cast<int>(y))) ||
         !sample.valid || sample.clipped ||
         inside_source_mask(x, y, sources, options.source_mask_radius)) {
        continue;
      }
      const int channel = channel_index(cfa.pattern().at(x, y));
      numerator[channel][y * width + x] = sample.value;
      denominator[channel][y * width + x] = 1.0;
    }
  }
  astrocfa::RgbImage diffuse(width, height);
  for(int channel = 0; channel < 3; ++channel) {
    numerator[channel] = gaussian_blur(numerator[channel], width, height,
                                       options.diffuse_sigma);
    denominator[channel] = gaussian_blur(denominator[channel], width, height,
                                         options.diffuse_sigma);
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      astrocfa::RgbPixel pixel;
      for(int channel = 0; channel < 3; ++channel) {
        const std::size_t index = y * width + x;
        const double value = denominator[channel][index] > 1.0e-8
                                 ? numerator[channel][index] /
                                       denominator[channel][index]
                                 : channel_value(fallback.pixel(x, y), channel);
        set_channel(pixel, channel, value);
      }
      diffuse.set_pixel(x, y, pixel);
    }
  }
  return diffuse;
}

double gaussian_atom(double x, double y, const SourceCandidate &source,
                     double sigma) {
  const double dx = x - source.x;
  const double dy = y - source.y;
  return std::exp(-0.5 * (dx * dx + dy * dy) / (sigma * sigma));
}

double profile_atom(double x, double y, const SourceCandidate &source,
                    const RadialProfile &profile) {
  const double dx = x - source.x;
  const double dy = y - source.y;
  const double coordinate = std::sqrt(dx * dx + dy * dy) / profile.bin_size;
  const std::size_t lower = static_cast<std::size_t>(std::floor(coordinate));
  if(lower >= profile.values.size()) {
    return 0.0;
  }
  const std::size_t upper = std::min(lower + 1U, profile.values.size() - 1U);
  const double fraction = coordinate - static_cast<double>(lower);
  return (1.0 - fraction) * profile.values[lower] +
         fraction * profile.values[upper];
}

double epsf_atom(double x, double y, const SourceCandidate &source,
                 const EpsfProfile &profile) {
  const double grid_x = (x - source.x + profile.radius) / profile.spacing;
  const double grid_y = (y - source.y + profile.radius) / profile.spacing;
  if(grid_x < 0.0 || grid_y < 0.0 ||
     grid_x > static_cast<double>(profile.size - 1U) ||
     grid_y > static_cast<double>(profile.size - 1U)) {
    return 0.0;
  }
  const std::size_t x0 = static_cast<std::size_t>(std::floor(grid_x));
  const std::size_t y0 = static_cast<std::size_t>(std::floor(grid_y));
  const std::size_t x1 = std::min(x0 + 1U, profile.size - 1U);
  const std::size_t y1 = std::min(y0 + 1U, profile.size - 1U);
  const double fx = grid_x - static_cast<double>(x0);
  const double fy = grid_y - static_cast<double>(y0);
  const double top = (1.0 - fx) * profile.values[y0 * profile.size + x0] +
                     fx * profile.values[y0 * profile.size + x1];
  const double bottom = (1.0 - fx) * profile.values[y1 * profile.size + x0] +
                        fx * profile.values[y1 * profile.size + x1];
  return (1.0 - fy) * top + fy * bottom;
}

double transformed_epsf_atom(double x, double y,
                             const SourceCandidate &source,
                             const EpsfProfile &profile, int channel,
                             const ChromaticEpsfTransform *transform) {
  if(transform == nullptr) {
    return epsf_atom(x, y, source, profile);
  }
  const double scale = transform->scale[channel];
  const double warped_x = source.x +
                          (x - source.x - transform->shift_x[channel]) / scale;
  const double warped_y = source.y +
                          (y - source.y - transform->shift_y[channel]) / scale;
  return epsf_atom(warped_x, warped_y, source, profile) / (scale * scale);
}

bool validation_sample(int x, int y) {
  return (x * 17 + y * 31) % 5 == 0;
}

SourceCandidate fit_source(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const SourceCandidate &detected, const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only = false) {
  SourceCandidate best = detected;
  double best_score = std::numeric_limits<double>::infinity();
  double null_score = 0.0;
  const int sample_radius =
      std::max(4, static_cast<int>(std::ceil(3.0 * options.maximum_source_sigma)));
  const int min_x = std::max(0, static_cast<int>(std::floor(detected.x)) - sample_radius);
  const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                             static_cast<int>(std::ceil(detected.x)) + sample_radius);
  const int min_y = std::max(0, static_cast<int>(std::floor(detected.y)) - sample_radius);
  const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                             static_cast<int>(std::ceil(detected.y)) + sample_radius);

  for(int y = min_y; y <= max_y; ++y) {
    for(int x = min_x; x <= max_x; ++x) {
      if(training_only && validation_sample(x, y)) {
        continue;
      }
      const astrocfa::CfaSample sample = cfa.sample_info(
          static_cast<std::size_t>(x), static_cast<std::size_t>(y));
      if(!sample.valid || sample.clipped) {
        continue;
      }
      const int channel = channel_index(cfa.pattern().at(
          static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
      const double residual = sample.value - channel_value(
          diffuse.pixel(static_cast<std::size_t>(x),
                        static_cast<std::size_t>(y)), channel);
      null_score += astrocfa::estimate_noise(sample.value, noise_model).weight *
                    residual * residual;
    }
  }

  for(double dy = -options.position_search_radius;
      dy <= options.position_search_radius + 1.0e-9;
      dy += options.position_search_step) {
    for(double dx = -options.position_search_radius;
        dx <= options.position_search_radius + 1.0e-9;
        dx += options.position_search_step) {
      SourceCandidate trial = detected;
      trial.x += dx;
      trial.y += dy;
      for(double sigma = options.minimum_source_sigma;
          sigma <= options.maximum_source_sigma + 1.0e-9;
          sigma += options.source_sigma_step) {
        std::array<double, 3> rhs{};
        std::array<double, 3> normal{};
        std::array<double, 3> lower_bound{};
        trial.sigma = sigma;
        for(int y = min_y; y <= max_y; ++y) {
          for(int x = min_x; x <= max_x; ++x) {
            if(training_only && validation_sample(x, y)) {
              continue;
            }
            const astrocfa::CfaSample sample = cfa.sample_info(
                static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            if(!sample.valid) {
              continue;
            }
            const int channel = channel_index(cfa.pattern().at(
                static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
            const double atom = gaussian_atom(x, y, trial, sigma);
            if(atom < 1.0e-6) {
              continue;
            }
            const double base = channel_value(
                diffuse.pixel(static_cast<std::size_t>(x),
                              static_cast<std::size_t>(y)),
                channel);
            if(sample.clipped) {
              lower_bound[channel] = std::max(
                  lower_bound[channel],
                  std::max(0.0, (static_cast<double>(sample.value) - base) / atom));
              continue;
            }
            const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight;
            rhs[channel] += weight * atom * (sample.value - base);
            normal[channel] += weight * atom * atom;
          }
        }
        for(int channel = 0; channel < 3; ++channel) {
          trial.amplitude[channel] = std::max(
              lower_bound[channel],
              normal[channel] > 0.0 ? std::max(0.0, rhs[channel] / normal[channel])
                                    : 0.0);
        }

        double score = 0.0;
        std::size_t samples = 0;
        for(int y = min_y; y <= max_y; ++y) {
          for(int x = min_x; x <= max_x; ++x) {
            if(training_only && validation_sample(x, y)) {
              continue;
            }
            const astrocfa::CfaSample sample = cfa.sample_info(
                static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            if(!sample.valid || sample.clipped) {
              continue;
            }
            const int channel = channel_index(cfa.pattern().at(
                static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
            const double prediction =
                channel_value(diffuse.pixel(static_cast<std::size_t>(x),
                                            static_cast<std::size_t>(y)),
                              channel) +
                trial.amplitude[channel] * gaussian_atom(x, y, trial, sigma);
            const double residual = sample.value - prediction;
            score += astrocfa::estimate_noise(sample.value, noise_model).weight *
                     residual * residual;
            samples += 1;
          }
        }
        if(samples > 0 && score < best_score) {
          best_score = score;
          best = trial;
        }
      }
    }
  }
  best.evidence = std::max(0.0, null_score - best_score);
  return best;
}

double source_validation_evidence(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const SourceCandidate &source, const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options) {
  const int radius =
      std::max(4, static_cast<int>(std::ceil(3.0 * options.maximum_source_sigma)));
  const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
  const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                             static_cast<int>(std::ceil(source.x)) + radius);
  const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
  const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                             static_cast<int>(std::ceil(source.y)) + radius);
  double evidence = 0.0;
  for(int y = min_y; y <= max_y; ++y) {
    for(int x = min_x; x <= max_x; ++x) {
      if(!validation_sample(x, y)) {
        continue;
      }
      const astrocfa::CfaSample sample = cfa.sample_info(
          static_cast<std::size_t>(x), static_cast<std::size_t>(y));
      if(!sample.valid || sample.clipped) {
        continue;
      }
      const int channel = channel_index(cfa.pattern().at(
          static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
      const double base = channel_value(
          diffuse.pixel(static_cast<std::size_t>(x),
                        static_cast<std::size_t>(y)), channel);
      const double prediction = base + source.amplitude[channel] *
          gaussian_atom(x, y, source, source.sigma);
      const double null_residual = sample.value - base;
      const double model_residual = sample.value - prediction;
      const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight;
      evidence += weight *
                  (null_residual * null_residual - model_residual * model_residual);
    }
  }
  return evidence;
}

std::vector<SourceCandidate> fit_shared_psf(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &initial,
    const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only = false) {
  std::vector<SourceCandidate> best = initial;
  double best_evidence = -1.0;
  for(double sigma = options.minimum_source_sigma;
      sigma <= options.maximum_source_sigma + 1.0e-9;
      sigma += options.source_sigma_step) {
    astrocfa::MorphologicalReconstructionOptions fixed = options;
    fixed.minimum_source_sigma = sigma;
    fixed.maximum_source_sigma = sigma;
    fixed.position_search_radius = std::min(0.25, options.position_search_radius);
    std::vector<SourceCandidate> fitted;
    fitted.reserve(initial.size());
    double evidence = 0.0;
    for(const SourceCandidate &source : initial) {
      fitted.push_back(
          fit_source(cfa, diffuse, source, noise_model, fixed, training_only));
      evidence += fitted.back().evidence;
    }
    if(evidence > best_evidence) {
      best_evidence = evidence;
      best = std::move(fitted);
    }
  }
  return best;
}

std::vector<SourceCandidate> fit_profile_amplitudes(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    std::vector<SourceCandidate> sources, const RadialProfile &profile,
    const astrocfa::NoiseModel &noise_model, bool training_only) {
  const int radius = static_cast<int>(
      std::ceil(profile.bin_size * static_cast<double>(profile.values.size() - 1U)));
  for(SourceCandidate &source : sources) {
    std::array<double, 3> rhs{};
    std::array<double, 3> normal{};
    const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
    const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                               static_cast<int>(std::ceil(source.x)) + radius);
    const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
    const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                               static_cast<int>(std::ceil(source.y)) + radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        if(training_only && validation_sample(x, y)) {
          continue;
        }
        const astrocfa::CfaSample sample = cfa.sample_info(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y));
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const int channel = channel_index(cfa.pattern().at(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
        const double atom = profile_atom(x, y, source, profile);
        const double base = channel_value(
            diffuse.pixel(static_cast<std::size_t>(x),
                          static_cast<std::size_t>(y)), channel);
        const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight;
        rhs[channel] += weight * atom * (sample.value - base);
        normal[channel] += weight * atom * atom;
      }
    }
    for(int channel = 0; channel < 3; ++channel) {
      source.amplitude[channel] = normal[channel] > 0.0
                                      ? std::max(0.0, rhs[channel] / normal[channel])
                                      : 0.0;
    }
  }
  return sources;
}

double log_determinant_spd(const std::vector<double> &matrix,
                           std::size_t dimension) {
  std::vector<double> lower(dimension * dimension, 0.0);
  double log_determinant = 0.0;
  for(std::size_t row = 0; row < dimension; ++row) {
    for(std::size_t column = 0; column <= row; ++column) {
      double value = matrix[row * dimension + column];
      for(std::size_t k = 0; k < column; ++k) {
        value -= lower[row * dimension + k] * lower[column * dimension + k];
      }
      if(row == column) {
        if(value <= 1.0e-12) {
          return -std::numeric_limits<double>::infinity();
        }
        lower[row * dimension + column] = std::sqrt(value);
        log_determinant += std::log(value);
      } else {
        lower[row * dimension + column] =
            value / lower[column * dimension + column];
      }
    }
  }
  return log_determinant;
}

ObservableSourceSelection select_observable_sources(
    const astrocfa::CfaFrame &cfa,
    const std::vector<SourceCandidate> &sources,
    const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only) {
  const std::size_t bins = static_cast<std::size_t>(
                               std::ceil(options.profile_radius /
                                         options.profile_bin_size)) +
                           1U;
  std::vector<std::vector<double>> information(
      sources.size(), std::vector<double>(bins * bins, 0.0));
  const int radius = static_cast<int>(std::ceil(options.profile_radius));
  for(std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
    const SourceCandidate &source = sources[source_index];
    const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
    const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                               static_cast<int>(std::ceil(source.x)) + radius);
    const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
    const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                               static_cast<int>(std::ceil(source.y)) + radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        if(training_only && validation_sample(x, y)) {
          continue;
        }
        const astrocfa::CfaSample sample = cfa.sample_info(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y));
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const int channel = channel_index(cfa.pattern().at(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
        const double amplitude = source.amplitude[channel];
        const double dx = static_cast<double>(x) - source.x;
        const double dy = static_cast<double>(y) - source.y;
        const double coordinate =
            std::sqrt(dx * dx + dy * dy) / options.profile_bin_size;
        const std::size_t lower = static_cast<std::size_t>(std::floor(coordinate));
        if(lower >= bins) {
          continue;
        }
        const std::size_t upper = std::min(lower + 1U, bins - 1U);
        const double fraction = coordinate - static_cast<double>(lower);
        const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight *
                              amplitude * amplitude;
        const double lower_basis = 1.0 - fraction;
        const double upper_basis = fraction;
        information[source_index][lower * bins + lower] +=
            weight * lower_basis * lower_basis;
        information[source_index][upper * bins + upper] +=
            weight * upper_basis * upper_basis;
        const double cross = weight * lower_basis * upper_basis;
        information[source_index][lower * bins + upper] += cross;
        information[source_index][upper * bins + lower] += cross;
      }
    }
  }

  std::vector<double> accumulated(bins * bins, 0.0);
  for(std::size_t bin = 0; bin < bins; ++bin) {
    accumulated[bin * bins + bin] = 1.0e-3;
  }
  const double prior_log_determinant = log_determinant_spd(accumulated, bins);
  double accumulated_log_determinant = prior_log_determinant;
  std::vector<bool> used(sources.size(), false);
  std::vector<SourceCandidate> selected;
  selected.reserve(std::min(sources.size(), options.maximum_profile_sources));
  while(selected.size() < std::min(sources.size(), options.maximum_profile_sources)) {
    std::size_t best_index = sources.size();
    double best_gain = -1.0;
    for(std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
      if(used[source_index]) {
        continue;
      }
      std::vector<double> candidate = accumulated;
      for(std::size_t element = 0; element < candidate.size(); ++element) {
        candidate[element] += information[source_index][element];
      }
      const double gain =
          log_determinant_spd(candidate, bins) - accumulated_log_determinant;
      if(gain > best_gain) {
        best_gain = gain;
        best_index = source_index;
      }
    }
    if(best_index == sources.size() || best_gain <= 0.0) {
      break;
    }
    used[best_index] = true;
    selected.push_back(sources[best_index]);
    for(std::size_t element = 0; element < accumulated.size(); ++element) {
      accumulated[element] += information[best_index][element];
    }
    accumulated_log_determinant = log_determinant_spd(accumulated, bins);
  }
  if(selected.empty()) {
    return {.sources = sources, .information_gain = 0.0};
  }
  return {.sources = std::move(selected),
          .information_gain = accumulated_log_determinant - prior_log_determinant};
}

RadialProfile learn_shared_profile(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    std::vector<SourceCandidate> &sources,
    const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only) {
  RadialProfile profile;
  profile.bin_size = options.profile_bin_size;
  const std::size_t bins = static_cast<std::size_t>(
                               std::ceil(options.profile_radius /
                                         options.profile_bin_size)) +
                           1U;
  profile.values.resize(bins);
  const double initial_sigma = sources.empty() ? 0.7 : sources.front().sigma;
  for(std::size_t bin = 0; bin < bins; ++bin) {
    const double radius = static_cast<double>(bin) * profile.bin_size;
    profile.values[bin] =
        std::exp(-0.5 * radius * radius / (initial_sigma * initial_sigma));
  }

  for(std::size_t iteration = 0; iteration < options.profile_iterations; ++iteration) {
    sources = fit_profile_amplitudes(cfa, diffuse, std::move(sources), profile,
                                     noise_model, training_only);
    std::vector<double> numerator(bins, 0.0);
    std::vector<double> denominator(bins, 0.0);
    const int radius = static_cast<int>(std::ceil(options.profile_radius));
    for(const SourceCandidate &source : sources) {
      const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
      const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                                 static_cast<int>(std::ceil(source.x)) + radius);
      const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
      const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                                 static_cast<int>(std::ceil(source.y)) + radius);
      for(int y = min_y; y <= max_y; ++y) {
        for(int x = min_x; x <= max_x; ++x) {
          if(training_only && validation_sample(x, y)) {
            continue;
          }
          const astrocfa::CfaSample sample = cfa.sample_info(
              static_cast<std::size_t>(x), static_cast<std::size_t>(y));
          if(!sample.valid || sample.clipped) {
            continue;
          }
          const int channel = channel_index(cfa.pattern().at(
              static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
          const double amplitude = source.amplitude[channel];
          if(amplitude <= 1.0e-5) {
            continue;
          }
          const double dx = static_cast<double>(x) - source.x;
          const double dy = static_cast<double>(y) - source.y;
          const std::size_t bin = static_cast<std::size_t>(
              std::llround(std::sqrt(dx * dx + dy * dy) / profile.bin_size));
          if(bin >= bins) {
            continue;
          }
          const double base = channel_value(
              diffuse.pixel(static_cast<std::size_t>(x),
                            static_cast<std::size_t>(y)), channel);
          const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight;
          numerator[bin] += weight * amplitude * (sample.value - base);
          denominator[bin] += weight * amplitude * amplitude;
        }
      }
    }
    std::vector<double> updated = profile.values;
    for(std::size_t bin = 0; bin < bins; ++bin) {
      if(denominator[bin] > 0.0) {
        updated[bin] = std::max(0.0, numerator[bin] / denominator[bin]);
      }
    }
    for(int pass = 0; pass < 2; ++pass) {
      std::vector<double> smooth = updated;
      for(std::size_t bin = 1; bin + 1U < bins; ++bin) {
        smooth[bin] = 0.25 * updated[bin - 1U] + 0.5 * updated[bin] +
                      0.25 * updated[bin + 1U];
      }
      updated = std::move(smooth);
    }
    const double normalization =
        std::max(1.0e-6, *std::max_element(updated.begin(), updated.end()));
    for(double &value : updated) {
      value = std::clamp(value / normalization, 0.0, 1.0);
    }
    for(std::size_t bin = 1; bin < bins; ++bin) {
      updated[bin] = std::min(updated[bin], updated[bin - 1U]);
    }
    profile.values = std::move(updated);
  }
  sources = fit_profile_amplitudes(cfa, diffuse, std::move(sources), profile,
                                   noise_model, training_only);
  return profile;
}

std::vector<SourceCandidate> fit_epsf_amplitudes(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    std::vector<SourceCandidate> sources, const EpsfProfile &profile,
    const astrocfa::NoiseModel &noise_model, bool training_only,
    const ChromaticEpsfTransform *transform = nullptr) {
  const int radius = static_cast<int>(std::ceil(profile.radius));
  for(SourceCandidate &source : sources) {
    std::array<double, 3> rhs{};
    std::array<double, 3> normal{};
    double null_score = 0.0;
    const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
    const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                               static_cast<int>(std::ceil(source.x)) + radius);
    const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
    const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                               static_cast<int>(std::ceil(source.y)) + radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        if(training_only && validation_sample(x, y)) {
          continue;
        }
        const astrocfa::CfaSample sample = cfa.sample_info(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y));
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const int channel = channel_index(cfa.pattern().at(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
        const double atom = transformed_epsf_atom(
            x, y, source, profile, channel, transform);
        const double base = channel_value(
            diffuse.pixel(static_cast<std::size_t>(x),
                          static_cast<std::size_t>(y)), channel);
        const double weight = astrocfa::estimate_noise(sample.value, noise_model).weight;
        rhs[channel] += weight * atom * (sample.value - base);
        normal[channel] += weight * atom * atom;
        const double residual = sample.value - base;
        null_score += weight * residual * residual;
      }
    }
    for(int channel = 0; channel < 3; ++channel) {
      source.amplitude[channel] = normal[channel] > 0.0
                                      ? std::max(0.0, rhs[channel] / normal[channel])
                                      : 0.0;
    }
    double model_score = 0.0;
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        if(training_only && validation_sample(x, y)) {
          continue;
        }
        const astrocfa::CfaSample sample = cfa.sample_info(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y));
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const int channel = channel_index(cfa.pattern().at(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
        const double base = channel_value(
            diffuse.pixel(static_cast<std::size_t>(x),
                          static_cast<std::size_t>(y)), channel);
        const double residual = sample.value - base -
                                source.amplitude[channel] *
                                    transformed_epsf_atom(
                                        x, y, source, profile, channel, transform);
        model_score += astrocfa::estimate_noise(sample.value, noise_model).weight *
                       residual * residual;
      }
    }
    source.evidence = std::max(0.0, null_score - model_score);
  }
  return sources;
}

std::vector<SourceCandidate> refine_epsf_sources(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &sources, const EpsfProfile &profile,
    const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only,
    const ChromaticEpsfTransform *transform = nullptr) {
  const double search_radius = std::min(0.375, options.position_search_radius);
  const double search_step = std::min(0.0625, options.position_search_step);
  std::vector<SourceCandidate> refined;
  refined.reserve(sources.size());
  for(const SourceCandidate &source : sources) {
    SourceCandidate best = source;
    double best_evidence = -1.0;
    for(double dy = -search_radius; dy <= search_radius + 1.0e-9;
        dy += search_step) {
      for(double dx = -search_radius; dx <= search_radius + 1.0e-9;
          dx += search_step) {
        SourceCandidate trial = source;
        trial.x += dx;
        trial.y += dy;
        std::vector<SourceCandidate> fitted = fit_epsf_amplitudes(
            cfa, diffuse, {trial}, profile, noise_model, training_only, transform);
        if(fitted.front().evidence > best_evidence) {
          best_evidence = fitted.front().evidence;
          best = fitted.front();
        }
      }
    }
    refined.push_back(best);
  }
  return refined;
}

double epsf_channel_score(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &sources, const EpsfProfile &profile,
    const astrocfa::NoiseModel &noise_model,
    const ChromaticEpsfTransform &transform, int target_channel,
    bool training_only) {
  double score = 0.0;
  const int radius = static_cast<int>(std::ceil(profile.radius));
  for(const SourceCandidate &source : sources) {
    const int min_x = std::max(0, static_cast<int>(std::floor(source.x)) - radius);
    const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                               static_cast<int>(std::ceil(source.x)) + radius);
    const int min_y = std::max(0, static_cast<int>(std::floor(source.y)) - radius);
    const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                               static_cast<int>(std::ceil(source.y)) + radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        if(training_only && validation_sample(x, y)) {
          continue;
        }
        const astrocfa::CfaSample sample = cfa.sample_info(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y));
        const int channel = channel_index(cfa.pattern().at(
            static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
        if(!sample.valid || sample.clipped || channel != target_channel) {
          continue;
        }
        const double prediction = channel_value(
                                      diffuse.pixel(static_cast<std::size_t>(x),
                                                    static_cast<std::size_t>(y)),
                                      channel) +
                                  source.amplitude[channel] *
                                      transformed_epsf_atom(
                                          x, y, source, profile, channel, &transform);
        const double residual = sample.value - prediction;
        score += astrocfa::estimate_noise(sample.value, noise_model).weight *
                 residual * residual;
      }
    }
  }
  return score;
}

ChromaticEpsfTransform estimate_chromatic_transform(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &sources, const EpsfProfile &profile,
    const astrocfa::NoiseModel &noise_model, bool training_only) {
  ChromaticEpsfTransform transform;
  constexpr std::array<double, 5> scales = {0.90, 0.95, 1.0, 1.05, 1.10};
  for(const int channel : {0, 2}) {
    double best_score = std::numeric_limits<double>::infinity();
    ChromaticEpsfTransform best = transform;
    for(const double scale : scales) {
      for(double shift_y = -0.25; shift_y <= 0.25 + 1.0e-9;
          shift_y += 0.0625) {
        for(double shift_x = -0.25; shift_x <= 0.25 + 1.0e-9;
            shift_x += 0.0625) {
          ChromaticEpsfTransform trial = transform;
          trial.scale[channel] = scale;
          trial.shift_x[channel] = shift_x;
          trial.shift_y[channel] = shift_y;
          const std::vector<SourceCandidate> fitted = fit_epsf_amplitudes(
              cfa, diffuse, sources, profile, noise_model, training_only, &trial);
          const double score = epsf_channel_score(
              cfa, diffuse, fitted, profile, noise_model, trial, channel,
              training_only);
          if(score < best_score) {
            best_score = score;
            best = trial;
          }
        }
      }
    }
    transform = best;
  }
  return transform;
}

astrocfa::RgbImage estimate_diffuse_after_epsf(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &fallback,
    const std::vector<SourceCandidate> &sources, const EpsfProfile &profile,
    const astrocfa::MorphologicalReconstructionOptions &options,
    const ChromaticEpsfTransform *transform = nullptr) {
  const std::size_t width = cfa.width();
  const std::size_t height = cfa.height();
  std::array<std::vector<double>, 3> numerator;
  std::array<std::vector<double>, 3> denominator;
  for(int channel = 0; channel < 3; ++channel) {
    numerator[channel].assign(width * height, 0.0);
    denominator[channel].assign(width * height, 0.0);
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const astrocfa::CfaSample sample = cfa.sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }
      const int channel = channel_index(cfa.pattern().at(x, y));
      numerator[channel][y * width + x] = sample.value;
      denominator[channel][y * width + x] = 1.0;
    }
  }
  double support = profile.radius;
  if(transform != nullptr) {
    for(int channel = 0; channel < 3; ++channel) {
      support = std::max(
          support, profile.radius * transform->scale[channel] +
                       std::hypot(transform->shift_x[channel],
                                  transform->shift_y[channel]));
    }
  }
  const int sample_radius = static_cast<int>(std::ceil(support)) + 1;
  for(const SourceCandidate &source : sources) {
    const int min_x =
        std::max(0, static_cast<int>(std::floor(source.x)) - sample_radius);
    const int max_x = std::min(static_cast<int>(width) - 1,
                               static_cast<int>(std::ceil(source.x)) + sample_radius);
    const int min_y =
        std::max(0, static_cast<int>(std::floor(source.y)) - sample_radius);
    const int max_y = std::min(static_cast<int>(height) - 1,
                               static_cast<int>(std::ceil(source.y)) + sample_radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        const std::size_t sx = static_cast<std::size_t>(x);
        const std::size_t sy = static_cast<std::size_t>(y);
        const astrocfa::CfaSample sample = cfa.sample_info(sx, sy);
        if(!sample.valid || sample.clipped) {
          continue;
        }
        const int channel = channel_index(cfa.pattern().at(sx, sy));
        numerator[channel][sy * width + sx] -=
            source.amplitude[channel] * transformed_epsf_atom(
                                            x, y, source, profile, channel,
                                            transform);
      }
    }
  }
  astrocfa::RgbImage diffuse(width, height);
  for(int channel = 0; channel < 3; ++channel) {
    for(double &value : numerator[channel]) {
      value = std::max(0.0, value);
    }
    numerator[channel] = gaussian_blur(numerator[channel], width, height,
                                       options.diffuse_sigma);
    denominator[channel] = gaussian_blur(denominator[channel], width, height,
                                         options.diffuse_sigma);
  }
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      astrocfa::RgbPixel pixel;
      for(int channel = 0; channel < 3; ++channel) {
        const std::size_t index = y * width + x;
        const double value = denominator[channel][index] > 1.0e-8
                                 ? numerator[channel][index] /
                                       denominator[channel][index]
                                 : channel_value(fallback.pixel(x, y), channel);
        set_channel(pixel, channel, value);
      }
      diffuse.set_pixel(x, y, pixel);
    }
  }
  return diffuse;
}

void normalize_epsf(EpsfProfile &profile) {
  double integral = 0.0;
  for(double value : profile.values) {
    integral += value;
  }
  integral *= profile.spacing * profile.spacing;
  if(integral <= 1.0e-12) {
    throw std::runtime_error("Cannot normalize an empty ePSF");
  }
  for(double &value : profile.values) {
    value /= integral;
  }
}

EpsfProfile learn_shared_epsf(
    const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &diffuse,
    std::vector<SourceCandidate> &sources,
    const astrocfa::NoiseModel &noise_model,
    const astrocfa::MorphologicalReconstructionOptions &options,
    bool training_only) {
  EpsfProfile profile;
  profile.spacing = 1.0 / static_cast<double>(options.epsf_oversampling);
  profile.radius = options.profile_radius;
  profile.size = 2U * static_cast<std::size_t>(
                          std::ceil(profile.radius / profile.spacing)) +
                 1U;
  profile.values.resize(profile.size * profile.size);
  const double initial_sigma = sources.empty() ? 0.7 : sources.front().sigma;
  for(std::size_t y = 0; y < profile.size; ++y) {
    for(std::size_t x = 0; x < profile.size; ++x) {
      const double dx = static_cast<double>(x) * profile.spacing - profile.radius;
      const double dy = static_cast<double>(y) * profile.spacing - profile.radius;
      profile.values[y * profile.size + x] =
          std::exp(-0.5 * (dx * dx + dy * dy) /
                   (initial_sigma * initial_sigma));
    }
  }
  normalize_epsf(profile);

  const int sample_radius = static_cast<int>(std::ceil(profile.radius));
  for(std::size_t iteration = 0; iteration < options.epsf_iterations; ++iteration) {
    sources = fit_epsf_amplitudes(cfa, diffuse, std::move(sources), profile,
                                  noise_model, training_only);
    std::vector<double> gradient(profile.values.size(), 0.0);
    std::vector<double> curvature(profile.values.size(), 0.0);
    for(const SourceCandidate &source : sources) {
      const int min_x =
          std::max(0, static_cast<int>(std::floor(source.x)) - sample_radius);
      const int max_x = std::min(static_cast<int>(cfa.width()) - 1,
                                 static_cast<int>(std::ceil(source.x)) + sample_radius);
      const int min_y =
          std::max(0, static_cast<int>(std::floor(source.y)) - sample_radius);
      const int max_y = std::min(static_cast<int>(cfa.height()) - 1,
                                 static_cast<int>(std::ceil(source.y)) + sample_radius);
      for(int y = min_y; y <= max_y; ++y) {
        for(int x = min_x; x <= max_x; ++x) {
          if(training_only && validation_sample(x, y)) {
            continue;
          }
          const astrocfa::CfaSample sample = cfa.sample_info(
              static_cast<std::size_t>(x), static_cast<std::size_t>(y));
          if(!sample.valid || sample.clipped) {
            continue;
          }
          const int channel = channel_index(cfa.pattern().at(
              static_cast<std::size_t>(x), static_cast<std::size_t>(y)));
          const double amplitude = source.amplitude[channel];
          if(amplitude <= 1.0e-8) {
            continue;
          }
          const double grid_x =
              (static_cast<double>(x) - source.x + profile.radius) /
              profile.spacing;
          const double grid_y =
              (static_cast<double>(y) - source.y + profile.radius) /
              profile.spacing;
          if(grid_x < 0.0 || grid_y < 0.0 ||
             grid_x > static_cast<double>(profile.size - 1U) ||
             grid_y > static_cast<double>(profile.size - 1U)) {
            continue;
          }
          const std::size_t x0 = static_cast<std::size_t>(std::floor(grid_x));
          const std::size_t y0 = static_cast<std::size_t>(std::floor(grid_y));
          const std::size_t x1 = std::min(x0 + 1U, profile.size - 1U);
          const std::size_t y1 = std::min(y0 + 1U, profile.size - 1U);
          const double fx = grid_x - static_cast<double>(x0);
          const double fy = grid_y - static_cast<double>(y0);
          const std::array<std::size_t, 4> indices = {
              y0 * profile.size + x0, y0 * profile.size + x1,
              y1 * profile.size + x0, y1 * profile.size + x1};
          const std::array<double, 4> basis = {
              (1.0 - fx) * (1.0 - fy), fx * (1.0 - fy),
              (1.0 - fx) * fy, fx * fy};
          const double base = channel_value(
              diffuse.pixel(static_cast<std::size_t>(x),
                            static_cast<std::size_t>(y)), channel);
          const double residual = sample.value - base -
                                  amplitude * epsf_atom(x, y, source, profile);
          const double weight =
              astrocfa::estimate_noise(sample.value, noise_model).weight;
          for(std::size_t corner = 0; corner < indices.size(); ++corner) {
            const double derivative = amplitude * basis[corner];
            gradient[indices[corner]] += weight * derivative * residual;
            curvature[indices[corner]] += weight * derivative * derivative;
          }
        }
      }
    }

    std::vector<double> updated = profile.values;
    for(std::size_t index = 0; index < updated.size(); ++index) {
      if(curvature[index] > 1.0e-12) {
        updated[index] = std::max(
            0.0, updated[index] + 0.35 * gradient[index] / curvature[index]);
      }
    }
    std::vector<double> smooth = updated;
    for(std::size_t y = 1; y + 1U < profile.size; ++y) {
      for(std::size_t x = 1; x + 1U < profile.size; ++x) {
        const double dx = static_cast<double>(x) * profile.spacing - profile.radius;
        const double dy = static_cast<double>(y) * profile.spacing - profile.radius;
        const std::size_t index = y * profile.size + x;
        if(dx * dx + dy * dy > profile.radius * profile.radius) {
          smooth[index] = 0.0;
          continue;
        }
        smooth[index] = 0.6 * updated[index] +
                        0.1 * (updated[index - 1U] + updated[index + 1U] +
                               updated[index - profile.size] +
                               updated[index + profile.size]);
      }
    }
    profile.values = std::move(smooth);
    normalize_epsf(profile);
  }
  sources = fit_epsf_amplitudes(cfa, diffuse, std::move(sources), profile,
                                noise_model, training_only);
  return profile;
}

astrocfa::RgbImage compose(const astrocfa::RgbImage &diffuse,
                           const std::vector<SourceCandidate> &sources);

double validation_score(const astrocfa::CfaFrame &cfa,
                        const astrocfa::RgbImage &model,
                        const astrocfa::NoiseModel &noise_model) {
  double score = 0.0;
  std::size_t samples = 0;
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      if(!validation_sample(static_cast<int>(x), static_cast<int>(y))) {
        continue;
      }
      const astrocfa::CfaSample sample = cfa.sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }
      const int channel = channel_index(cfa.pattern().at(x, y));
      const double residual = sample.value - channel_value(model.pixel(x, y), channel);
      score += astrocfa::estimate_noise(sample.value, noise_model).weight *
               residual * residual;
      samples += 1;
    }
  }
  return samples > 0 ? score / static_cast<double>(samples)
                     : std::numeric_limits<double>::infinity();
}

std::size_t validation_sample_count(const astrocfa::CfaFrame &cfa) {
  std::size_t samples = 0;
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const astrocfa::CfaSample sample = cfa.sample_info(x, y);
      if(validation_sample(static_cast<int>(x), static_cast<int>(y)) &&
         sample.valid && !sample.clipped) {
        samples += 1U;
      }
    }
  }
  return samples;
}

astrocfa::RgbImage compose_profile(
    const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &sources,
    const RadialProfile &profile) {
  astrocfa::RgbImage output(diffuse.width(), diffuse.height());
  for(std::size_t y = 0; y < diffuse.height(); ++y) {
    for(std::size_t x = 0; x < diffuse.width(); ++x) {
      astrocfa::RgbPixel pixel = diffuse.pixel(x, y);
      for(const SourceCandidate &source : sources) {
        const double atom = profile_atom(x, y, source, profile);
        pixel.r += static_cast<float>(source.amplitude[0] * atom);
        pixel.g += static_cast<float>(source.amplitude[1] * atom);
        pixel.b += static_cast<float>(source.amplitude[2] * atom);
      }
      pixel.r = std::clamp(pixel.r, 0.0F, 1.25F);
      pixel.g = std::clamp(pixel.g, 0.0F, 1.25F);
      pixel.b = std::clamp(pixel.b, 0.0F, 1.25F);
      output.set_pixel(x, y, pixel);
    }
  }
  return output;
}

astrocfa::RgbImage compose_epsf(
    const astrocfa::RgbImage &diffuse,
    const std::vector<SourceCandidate> &sources,
    const EpsfProfile &profile,
    const ChromaticEpsfTransform *transform = nullptr) {
  astrocfa::RgbImage output = diffuse;
  double support = profile.radius;
  if(transform != nullptr) {
    for(int channel = 0; channel < 3; ++channel) {
      support = std::max(
          support, profile.radius * transform->scale[channel] +
                       std::hypot(transform->shift_x[channel],
                                  transform->shift_y[channel]));
    }
  }
  const int sample_radius = static_cast<int>(std::ceil(support)) + 1;
  for(const SourceCandidate &source : sources) {
    const int min_x =
        std::max(0, static_cast<int>(std::floor(source.x)) - sample_radius);
    const int max_x =
        std::min(static_cast<int>(diffuse.width()) - 1,
                 static_cast<int>(std::ceil(source.x)) + sample_radius);
    const int min_y =
        std::max(0, static_cast<int>(std::floor(source.y)) - sample_radius);
    const int max_y =
        std::min(static_cast<int>(diffuse.height()) - 1,
                 static_cast<int>(std::ceil(source.y)) + sample_radius);
    for(int y = min_y; y <= max_y; ++y) {
      for(int x = min_x; x <= max_x; ++x) {
        astrocfa::RgbPixel pixel = output.pixel(static_cast<std::size_t>(x),
                                                static_cast<std::size_t>(y));
        pixel.r += static_cast<float>(
            source.amplitude[0] *
            transformed_epsf_atom(x, y, source, profile, 0, transform));
        pixel.g += static_cast<float>(
            source.amplitude[1] *
            transformed_epsf_atom(x, y, source, profile, 1, transform));
        pixel.b += static_cast<float>(
            source.amplitude[2] *
            transformed_epsf_atom(x, y, source, profile, 2, transform));
        output.set_pixel(static_cast<std::size_t>(x), static_cast<std::size_t>(y),
                         pixel);
      }
    }
  }
  for(std::size_t y = 0; y < diffuse.height(); ++y) {
    for(std::size_t x = 0; x < diffuse.width(); ++x) {
      astrocfa::RgbPixel pixel = output.pixel(x, y);
      pixel.r = std::clamp(pixel.r, 0.0F, 1.25F);
      pixel.g = std::clamp(pixel.g, 0.0F, 1.25F);
      pixel.b = std::clamp(pixel.b, 0.0F, 1.25F);
      output.set_pixel(x, y, pixel);
    }
  }
  return output;
}

astrocfa::RgbImage compose(const astrocfa::RgbImage &diffuse,
                           const std::vector<SourceCandidate> &sources) {
  astrocfa::RgbImage output(diffuse.width(), diffuse.height());
  for(std::size_t y = 0; y < diffuse.height(); ++y) {
    for(std::size_t x = 0; x < diffuse.width(); ++x) {
      astrocfa::RgbPixel pixel = diffuse.pixel(x, y);
      for(const SourceCandidate &source : sources) {
        const double atom = gaussian_atom(x, y, source, source.sigma);
        pixel.r += static_cast<float>(source.amplitude[0] * atom);
        pixel.g += static_cast<float>(source.amplitude[1] * atom);
        pixel.b += static_cast<float>(source.amplitude[2] * atom);
      }
      pixel.r = std::clamp(pixel.r, 0.0F, 1.25F);
      pixel.g = std::clamp(pixel.g, 0.0F, 1.25F);
      pixel.b = std::clamp(pixel.b, 0.0F, 1.25F);
      output.set_pixel(x, y, pixel);
    }
  }
  return output;
}

} // namespace

namespace astrocfa {

const char *morphological_model_name(MorphologicalModel model) {
  switch(model) {
  case MorphologicalModel::independent_gaussian:
    return "independent-gaussian";
  case MorphologicalModel::shared_gaussian:
    return "shared-gaussian";
  case MorphologicalModel::shared_profile:
    return "shared-profile";
  case MorphologicalModel::shared_epsf:
    return "shared-epsf";
  case MorphologicalModel::shared_chromatic_epsf:
    return "shared-chromatic-epsf";
  }
  return "unknown";
}

MorphologicalReconstructionResult reconstruct_morphological_cfa_detailed(
    const CfaFrame &cfa, const NoiseModel &noise_model,
    MorphologicalReconstructionOptions options) {
  if(options.diffuse_sigma <= 0.0 || options.detection_sigma <= 0.0 ||
     options.source_mask_radius <= 0.0 || options.minimum_source_sigma <= 0.0 ||
     options.maximum_source_sigma < options.minimum_source_sigma ||
     options.source_sigma_step <= 0.0 || options.position_search_radius < 0.0 ||
     options.position_search_step <= 0.0 || options.maximum_sources == 0 ||
     options.profile_bin_size <= 0.0 || options.profile_radius <= 0.0 ||
     options.profile_iterations == 0 || options.maximum_profile_sources == 0 ||
     options.epsf_oversampling == 0 || options.epsf_oversampling > 16 ||
     options.epsf_iterations == 0) {
    throw std::invalid_argument("Invalid morphological CFA reconstruction options");
  }

  const DemosaicResult initial = reconstruct_frequency_guided(cfa, noise_model);
  std::vector<SourceCandidate> sources =
      detect_sources(cfa, noise_model, options);
  const std::vector<SourceCandidate> detected_sources = sources;
  MorphologicalReconstructionStats stats;
  stats.detected_sources = sources.size();
  RgbImage training_diffuse =
      estimate_diffuse(cfa, initial.image, sources, options, true);
  std::vector<SourceCandidate> validated_sources;
  validated_sources.reserve(sources.size());
  const double minimum_evidence =
      options.detection_sigma * options.detection_sigma;
  for(const SourceCandidate &source : sources) {
    const SourceCandidate training =
        fit_source(cfa, training_diffuse, source, noise_model, options, true);
    if(training.evidence >= minimum_evidence &&
       source_validation_evidence(cfa, training_diffuse, training, noise_model,
                                  options) > 0.0) {
      validated_sources.push_back(training);
    }
  }
  sources = std::move(validated_sources);
  stats.validated_sources = sources.size();

  MorphologicalModel selected_model = MorphologicalModel::independent_gaussian;
  if(options.share_psf_across_sources && sources.size() > 1U) {
    std::vector<SourceCandidate> independent_training = sources;
    for(SourceCandidate &source : independent_training) {
      source = fit_source(cfa, training_diffuse, source, noise_model, options, true);
    }
    const std::vector<SourceCandidate> shared_training = fit_shared_psf(
        cfa, training_diffuse, sources, noise_model, options, true);
    ObservableSourceSelection profile_selection =
        select_observable_sources(cfa, shared_training, noise_model, options, true);
    stats.profile_sources = profile_selection.sources.size();
    stats.profile_information_gain = profile_selection.information_gain;
    const RadialProfile profile_training = learn_shared_profile(
        cfa, training_diffuse, profile_selection.sources, noise_model, options, true);
    std::vector<SourceCandidate> profile_training_sources = fit_profile_amplitudes(
        cfa, training_diffuse, shared_training, profile_training, noise_model, true);
    stats.independent_validation = validation_score(
        cfa, compose(training_diffuse, independent_training), noise_model);
    stats.shared_validation = validation_score(
        cfa, compose(training_diffuse, shared_training), noise_model);
    stats.profile_validation = validation_score(
        cfa, compose_profile(training_diffuse, profile_training_sources,
                             profile_training), noise_model);
    if(options.enable_epsf) {
      std::vector<SourceCandidate> epsf_learning_sources = shared_training;
      const EpsfProfile epsf_training = learn_shared_epsf(
          cfa, training_diffuse, epsf_learning_sources, noise_model, options, true);
      std::vector<SourceCandidate> epsf_training_sources;
      epsf_training_sources.reserve(detected_sources.size());
      for(const SourceCandidate &candidate : detected_sources) {
        epsf_training_sources.push_back(
            fit_source(cfa, training_diffuse, candidate, noise_model, options, true));
      }
      epsf_training_sources = fit_epsf_amplitudes(
          cfa, training_diffuse, std::move(epsf_training_sources), epsf_training,
          noise_model, true);
      epsf_training_sources = refine_epsf_sources(
          cfa, training_diffuse, epsf_training_sources, epsf_training,
          noise_model, options, true);
      std::erase_if(epsf_training_sources, [&](const SourceCandidate &source) {
        return source.evidence < minimum_evidence;
      });
      stats.epsf_validation = validation_score(
          cfa, compose_epsf(training_diffuse, epsf_training_sources, epsf_training),
          noise_model);
      if(options.enable_chromatic_epsf) {
        const ChromaticEpsfTransform chromatic_training =
            estimate_chromatic_transform(cfa, training_diffuse,
                                         epsf_training_sources, epsf_training,
                                         noise_model, true);
        const std::vector<SourceCandidate> chromatic_training_sources =
            fit_epsf_amplitudes(cfa, training_diffuse, epsf_training_sources,
                                epsf_training, noise_model, true,
                                &chromatic_training);
        stats.chromatic_epsf_validation = validation_score(
            cfa, compose_epsf(training_diffuse, chromatic_training_sources,
                              epsf_training, &chromatic_training),
            noise_model);
      }
    }
    const std::size_t holdout_samples = validation_sample_count(cfa);
    stats.chromatic_complexity_penalty =
        holdout_samples > 1U
            ? 6.0 * std::log(static_cast<double>(holdout_samples)) /
                  static_cast<double>(holdout_samples)
            : std::numeric_limits<double>::infinity();
    const double penalized_chromatic_validation =
        stats.chromatic_epsf_validation + stats.chromatic_complexity_penalty;
    if(penalized_chromatic_validation < stats.independent_validation &&
       penalized_chromatic_validation < stats.shared_validation &&
       penalized_chromatic_validation < stats.profile_validation &&
       penalized_chromatic_validation < stats.epsf_validation) {
      selected_model = MorphologicalModel::shared_chromatic_epsf;
    } else if(stats.epsf_validation < stats.independent_validation &&
       stats.epsf_validation < stats.shared_validation &&
       stats.epsf_validation < stats.profile_validation) {
      selected_model = MorphologicalModel::shared_epsf;
    } else if(stats.profile_validation < stats.independent_validation &&
              stats.profile_validation < stats.shared_validation) {
      selected_model = MorphologicalModel::shared_profile;
    } else if(stats.shared_validation < stats.independent_validation) {
      selected_model = MorphologicalModel::shared_gaussian;
    }
  }
  stats.selected_model = selected_model;

  // Model selection is complete. Reintroduce held-out sensels only now and
  // refit the selected family against all valid measurements.
  RgbImage diffuse = estimate_diffuse(
      cfa, initial.image,
      selected_model == MorphologicalModel::shared_epsf ||
              selected_model == MorphologicalModel::shared_chromatic_epsf
          ? detected_sources
          : sources,
      options);
  RadialProfile selected_profile;
  EpsfProfile selected_epsf;
  ChromaticEpsfTransform selected_chromatic;
  if(selected_model == MorphologicalModel::shared_epsf ||
     selected_model == MorphologicalModel::shared_chromatic_epsf) {
    std::vector<SourceCandidate> epsf_learning_sources =
        fit_shared_psf(cfa, diffuse, sources, noise_model, options);
    selected_epsf = learn_shared_epsf(cfa, diffuse, epsf_learning_sources,
                                      noise_model, options, false);
    const ChromaticEpsfTransform *transform = nullptr;
    if(selected_model == MorphologicalModel::shared_chromatic_epsf) {
      selected_chromatic = estimate_chromatic_transform(
          cfa, diffuse, epsf_learning_sources, selected_epsf, noise_model, false);
      stats.red_epsf_shift_x = selected_chromatic.shift_x[0];
      stats.red_epsf_shift_y = selected_chromatic.shift_y[0];
      stats.red_epsf_scale = selected_chromatic.scale[0];
      stats.blue_epsf_shift_x = selected_chromatic.shift_x[2];
      stats.blue_epsf_shift_y = selected_chromatic.shift_y[2];
      stats.blue_epsf_scale = selected_chromatic.scale[2];
      transform = &selected_chromatic;
    }
    sources.clear();
    sources.reserve(detected_sources.size());
    for(const SourceCandidate &candidate : detected_sources) {
      sources.push_back(
          fit_source(cfa, diffuse, candidate, noise_model, options, false));
    }
    sources = fit_epsf_amplitudes(cfa, diffuse, std::move(sources), selected_epsf,
                                  noise_model, false, transform);
    sources = refine_epsf_sources(cfa, diffuse, sources, selected_epsf,
                                  noise_model, options, false, transform);
    for(int iteration = 0; iteration < 2; ++iteration) {
      std::erase_if(sources, [&](const SourceCandidate &source) {
        return source.evidence < minimum_evidence;
      });
      diffuse = estimate_diffuse_after_epsf(cfa, diffuse, sources, selected_epsf,
                                            options, transform);
      sources = refine_epsf_sources(cfa, diffuse, sources, selected_epsf,
                                    noise_model, options, false, transform);
    }
  } else if(selected_model == MorphologicalModel::shared_profile) {
    sources = fit_shared_psf(cfa, diffuse, sources, noise_model, options);
    ObservableSourceSelection profile_selection =
        select_observable_sources(cfa, sources, noise_model, options, false);
    selected_profile = learn_shared_profile(cfa, diffuse, profile_selection.sources,
                                            noise_model, options, false);
    sources = fit_profile_amplitudes(cfa, diffuse, std::move(sources),
                                     selected_profile, noise_model, false);
  } else if(selected_model == MorphologicalModel::shared_gaussian) {
    sources = fit_shared_psf(cfa, diffuse, sources, noise_model, options);
  } else {
    for(SourceCandidate &source : sources) {
      source = fit_source(cfa, diffuse, source, noise_model, options);
    }
  }
  std::erase_if(sources, [&](const SourceCandidate &source) {
    return source.evidence < minimum_evidence;
  });
  stats.reconstructed_sources = sources.size();
  RgbImage output = selected_model == MorphologicalModel::shared_epsf ||
                            selected_model == MorphologicalModel::shared_chromatic_epsf
                        ? compose_epsf(
                              diffuse, sources, selected_epsf,
                              selected_model == MorphologicalModel::shared_chromatic_epsf
                                  ? &selected_chromatic
                                  : nullptr)
                    : selected_model == MorphologicalModel::shared_profile
                        ? compose_profile(diffuse, sources, selected_profile)
                        : compose(diffuse, sources);
  return MorphologicalReconstructionResult{
      .reconstruction = DemosaicResult{
          .image = output,
          .residual = compute_remosaic_residual(cfa, output),
          .noise_weighted_residual =
              compute_noise_weighted_remosaic_residual(cfa, output, noise_model),
      },
      .stats = stats,
  };
}

DemosaicResult reconstruct_morphological_cfa(
    const CfaFrame &cfa, const NoiseModel &noise_model,
    MorphologicalReconstructionOptions options) {
  return reconstruct_morphological_cfa_detailed(cfa, noise_model, options)
      .reconstruction;
}

} // namespace astrocfa
