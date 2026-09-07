#include "astrocfa/calibration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int phase_index(astrocfa::CfaColor color) {
  return static_cast<int>(color);
}

void require_compatible(const astrocfa::CfaFrame &light,
                        const astrocfa::CfaFrame *master,
                        const char *name) {
  if(master == nullptr) {
    return;
  }
  if(master->width() != light.width() || master->height() != light.height()) {
    throw std::invalid_argument(std::string(name) + " dimensions do not match light");
  }
  for(std::size_t y = 0; y < 2 && y < light.height(); ++y) {
    for(std::size_t x = 0; x < 2 && x < light.width(); ++x) {
      if(master->pattern().at(x, y) != light.pattern().at(x, y)) {
        throw std::invalid_argument(std::string(name) + " Bayer phase does not match light");
      }
    }
  }
}

void compute_flat_phase_means(const astrocfa::CfaFrame &light,
                              const astrocfa::CfaFrame *flat,
                              double means[4]) {
  for(int i = 0; i < 4; ++i) {
    means[i] = 1.0;
  }
  if(flat == nullptr) {
    return;
  }

  double sums[4] = {0.0, 0.0, 0.0, 0.0};
  std::size_t counts[4] = {0, 0, 0, 0};
  for(std::size_t y = 0; y < light.height(); ++y) {
    for(std::size_t x = 0; x < light.width(); ++x) {
      const astrocfa::CfaSample sample = flat->sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }
      const int index = phase_index(light.pattern().at(x, y));
      sums[index] += sample.value;
      counts[index] += 1;
    }
  }

  for(int i = 0; i < 4; ++i) {
    if(counts[i] > 0 && sums[i] > 0.0) {
      means[i] = sums[i] / static_cast<double>(counts[i]);
    }
  }
}

void require_frame_compatible(const astrocfa::CfaFrame &reference,
                              const astrocfa::CfaFrame &frame,
                              const char *name) {
  if(frame.width() != reference.width() || frame.height() != reference.height()) {
    throw std::invalid_argument(std::string(name) + " dimensions do not match");
  }
  for(std::size_t y = 0; y < 2 && y < reference.height(); ++y) {
    for(std::size_t x = 0; x < 2 && x < reference.width(); ++x) {
      if(frame.pattern().at(x, y) != reference.pattern().at(x, y)) {
        throw std::invalid_argument(std::string(name) + " Bayer phase does not match");
      }
    }
  }
}

double median_in_place(std::vector<double> &values) {
  if(values.empty()) {
    return 0.0;
  }
  const std::size_t middle = values.size() / 2U;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle),
                   values.end());
  const double upper = values[middle];
  if(values.size() % 2U == 1U) {
    return upper;
  }
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle - 1U),
                   values.end());
  return 0.5 * (values[middle - 1U] + upper);
}

double robust_pixel_value(std::vector<double> values, const astrocfa::MasterBuildOptions &options,
                          std::size_t &rejected) {
  if(values.empty()) {
    return 0.0;
  }
  if(values.size() < options.min_clip_samples || options.sigma_clip <= 0.0) {
    return median_in_place(values);
  }

  std::vector<double> work = values;
  const double median = median_in_place(work);
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for(double value : values) {
    deviations.push_back(std::abs(value - median));
  }
  const double mad = median_in_place(deviations);
  const double sigma = 1.4826 * mad;
  if(sigma <= 1.0e-12) {
    return median;
  }

  std::vector<double> kept;
  kept.reserve(values.size());
  const double threshold = options.sigma_clip * sigma;
  for(double value : values) {
    if(std::abs(value - median) <= threshold) {
      kept.push_back(value);
    } else {
      rejected += 1;
    }
  }

  if(kept.empty()) {
    return median;
  }
  return median_in_place(kept);
}

} // namespace

namespace astrocfa {

CalibrationResult calibrate_cfa(const CfaFrame &light, CalibrationInputs inputs) {
  require_compatible(light, inputs.bias, "Bias master");
  require_compatible(light, inputs.dark, "Dark master");
  require_compatible(light, inputs.flat, "Flat master");

  CalibrationResult result{
      .cfa = CfaFrame(light.width(), light.height(), light.pattern()),
  };
  compute_flat_phase_means(light, inputs.flat, result.stats.flat_phase_mean);

  for(std::size_t y = 0; y < light.height(); ++y) {
    for(std::size_t x = 0; x < light.width(); ++x) {
      const CfaSample light_sample = light.sample_info(x, y);
      CfaSample corrected = light_sample;
      result.stats.samples += 1;
      result.stats.invalid_samples += light_sample.valid ? 0U : 1U;
      result.stats.clipped_samples += light_sample.clipped ? 1U : 0U;
      if(light_sample.valid && !light_sample.clipped) {
        result.stats.mean_before += light_sample.value;
      }

      double value = light_sample.value;
      if(light_sample.valid) {
        if(inputs.dark != nullptr) {
          const CfaSample dark_sample = inputs.dark->sample_info(x, y);
          if(dark_sample.valid && !dark_sample.clipped) {
            value -= dark_sample.value;
            result.stats.mean_dark_subtracted += dark_sample.value;
          } else {
            corrected.valid = false;
          }
        }
        if(inputs.bias != nullptr &&
           (inputs.dark == nullptr || !inputs.options.dark_includes_bias)) {
          const CfaSample bias_sample = inputs.bias->sample_info(x, y);
          if(bias_sample.valid && !bias_sample.clipped) {
            value -= bias_sample.value;
            result.stats.mean_bias_subtracted += bias_sample.value;
          } else {
            corrected.valid = false;
          }
        }
        if(inputs.flat != nullptr) {
          const CfaSample flat_sample = inputs.flat->sample_info(x, y);
          if(flat_sample.valid && !flat_sample.clipped) {
            const int index = phase_index(light.pattern().at(x, y));
            const double normalized_flat =
                result.stats.flat_phase_mean[index] > 0.0
                    ? flat_sample.value / result.stats.flat_phase_mean[index]
                    : 1.0;
            const double divisor = std::max(normalized_flat, inputs.options.flat_floor);
            result.stats.flat_floor_samples +=
                normalized_flat < inputs.options.flat_floor ? 1U : 0U;
            value /= divisor;
          } else {
            corrected.valid = false;
          }
        }
      }

      corrected.value = static_cast<float>(std::clamp(value, 0.0, 1.25));
      result.cfa.set_sample(x, y, corrected);
      if(corrected.valid && !corrected.clipped) {
        result.stats.mean_after += corrected.value;
      }
    }
  }

  const double usable = static_cast<double>(
      result.stats.samples - result.stats.invalid_samples - result.stats.clipped_samples);
  if(usable > 0.0) {
    result.stats.mean_before /= usable;
    result.stats.mean_after /= usable;
    result.stats.mean_bias_subtracted /= usable;
    result.stats.mean_dark_subtracted /= usable;
  }

  return result;
}

MasterBuildResult build_master_cfa(const std::vector<CfaFrame> &frames,
                                   MasterBuildOptions options) {
  if(frames.empty()) {
    throw std::invalid_argument("Cannot build a master from zero frames");
  }

  const CfaFrame &reference = frames.front();
  for(std::size_t i = 1; i < frames.size(); ++i) {
    require_frame_compatible(reference, frames[i], "Master input frame");
  }

  MasterBuildResult result{
      .cfa = CfaFrame(reference.width(), reference.height(), reference.pattern()),
  };
  result.stats.frames = frames.size();
  result.stats.samples = reference.width() * reference.height();

  std::vector<double> values;
  values.reserve(frames.size());
  for(std::size_t y = 0; y < reference.height(); ++y) {
    for(std::size_t x = 0; x < reference.width(); ++x) {
      values.clear();
      bool any_clipped = false;
      for(const CfaFrame &frame : frames) {
        const CfaSample sample = frame.sample_info(x, y);
        any_clipped = any_clipped || sample.clipped;
        if(sample.valid && !sample.clipped) {
          values.push_back(sample.value);
        }
      }

      CfaSample master_sample;
      master_sample.valid = !values.empty();
      master_sample.clipped = any_clipped && values.empty();
      if(master_sample.valid) {
        master_sample.value =
            static_cast<float>(robust_pixel_value(values, options,
                                                  result.stats.rejected_samples));
        result.stats.mean += master_sample.value;
      } else {
        result.stats.invalid_output_samples += 1;
      }
      result.cfa.set_sample(x, y, master_sample);
    }
  }

  const double valid_samples = static_cast<double>(
      result.stats.samples - result.stats.invalid_output_samples);
  if(valid_samples > 0.0) {
    result.stats.mean /= valid_samples;
  }

  return result;
}

} // namespace astrocfa
