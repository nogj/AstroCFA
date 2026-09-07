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
                              double means[4],
                              const astrocfa::DefectMap *defects = nullptr) {
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
      if(defects != nullptr && defects->defective(x, y)) {
        continue;
      }
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

void collect_same_phase_neighbors(const astrocfa::CfaFrame &frame,
                                  std::size_t x, std::size_t y,
                                  std::size_t radius, std::vector<double> &values,
                                  const astrocfa::DefectMap *defects = nullptr) {
  values.clear();
  const auto signed_x = static_cast<std::ptrdiff_t>(x);
  const auto signed_y = static_cast<std::ptrdiff_t>(y);
  const auto width = static_cast<std::ptrdiff_t>(frame.width());
  const auto height = static_cast<std::ptrdiff_t>(frame.height());
  const auto signed_radius = static_cast<std::ptrdiff_t>(radius);
  for(std::ptrdiff_t phase_y = -signed_radius; phase_y <= signed_radius; ++phase_y) {
    for(std::ptrdiff_t phase_x = -signed_radius; phase_x <= signed_radius; ++phase_x) {
      if(phase_x == 0 && phase_y == 0) {
        continue;
      }
      const std::ptrdiff_t neighbor_x = signed_x + phase_x * 2;
      const std::ptrdiff_t neighbor_y = signed_y + phase_y * 2;
      if(neighbor_x < 0 || neighbor_y < 0 || neighbor_x >= width ||
         neighbor_y >= height) {
        continue;
      }
      const auto nx = static_cast<std::size_t>(neighbor_x);
      const auto ny = static_cast<std::size_t>(neighbor_y);
      if(defects != nullptr && defects->defective(nx, ny)) {
        continue;
      }
      const astrocfa::CfaSample sample = frame.sample_info(nx, ny);
      if(sample.valid && !sample.clipped) {
        values.push_back(sample.value);
      }
    }
  }
}

astrocfa::DefectMap detect_sensor_defects(
    const astrocfa::CfaFrame &light, const astrocfa::CfaFrame *dark,
    const astrocfa::CfaFrame *flat,
    const astrocfa::CosmeticCorrectionOptions &options,
    astrocfa::CalibrationStats &stats) {
  astrocfa::DefectMap defects(light.width(), light.height());
  if(!options.enabled || (dark == nullptr && flat == nullptr)) {
    return defects;
  }
  if(options.detection_radius == 0 || options.repair_radius == 0 ||
     options.min_neighbors == 0) {
    throw std::invalid_argument("Cosmetic correction radii and neighbor count must be positive");
  }

  const std::size_t neighborhood_width = options.detection_radius * 2U + 1U;
  std::vector<double> neighbors;
  std::vector<double> work;
  std::vector<double> deviations;
  neighbors.reserve(neighborhood_width * neighborhood_width - 1U);
  work.reserve(neighbors.capacity());
  deviations.reserve(neighbors.capacity());
  for(std::size_t y = 0; y < light.height(); ++y) {
    for(std::size_t x = 0; x < light.width(); ++x) {
      if(dark != nullptr) {
        const astrocfa::CfaSample sample = dark->sample_info(x, y);
        if(!sample.valid) {
          defects.add(x, y, astrocfa::SensorDefect::invalid_master);
        } else if(sample.clipped) {
          defects.add(x, y, astrocfa::SensorDefect::hot);
        } else {
          collect_same_phase_neighbors(*dark, x, y, options.detection_radius,
                                       neighbors);
          if(neighbors.size() >= options.min_neighbors) {
            work.assign(neighbors.begin(), neighbors.end());
            const double local_median = median_in_place(work);
            deviations.clear();
            for(double value : neighbors) {
              deviations.push_back(std::abs(value - local_median));
            }
            const double robust_sigma = 1.4826 * median_in_place(deviations);
            const double threshold =
                std::max(options.hot_min_excess, options.hot_sigma * robust_sigma);
            if(static_cast<double>(sample.value) - local_median > threshold) {
              defects.add(x, y, astrocfa::SensorDefect::hot);
            }
          }
        }
      }

      if(flat != nullptr) {
        const astrocfa::CfaSample sample = flat->sample_info(x, y);
        if(!sample.valid) {
          defects.add(x, y, astrocfa::SensorDefect::invalid_master);
        } else if(!sample.clipped) {
          collect_same_phase_neighbors(*flat, x, y, options.detection_radius,
                                       neighbors);
          if(neighbors.size() >= options.min_neighbors) {
            const double local_median = median_in_place(neighbors);
            const double response = local_median > 0.0
                                        ? static_cast<double>(sample.value) / local_median
                                        : 1.0;
            if(response < options.dead_response_ratio &&
               local_median - static_cast<double>(sample.value) >
                   options.dead_min_deficit) {
              defects.add(x, y, astrocfa::SensorDefect::dead);
            }
          }
        }
      }
    }
  }

  for(std::size_t y = 0; y < light.height(); ++y) {
    for(std::size_t x = 0; x < light.width(); ++x) {
      stats.hot_pixels += defects.has(x, y, astrocfa::SensorDefect::hot) ? 1U : 0U;
      stats.dead_pixels += defects.has(x, y, astrocfa::SensorDefect::dead) ? 1U : 0U;
      stats.invalid_master_pixels +=
          defects.has(x, y, astrocfa::SensorDefect::invalid_master) ? 1U : 0U;
    }
  }
  return defects;
}

void repair_sensor_defects(astrocfa::CfaFrame &cfa,
                           const astrocfa::DefectMap &defects,
                           const astrocfa::CosmeticCorrectionOptions &options,
                           astrocfa::CalibrationStats &stats) {
  const astrocfa::CfaFrame measured = cfa;
  const std::size_t neighborhood_width = options.repair_radius * 2U + 1U;
  std::vector<double> neighbors;
  neighbors.reserve(neighborhood_width * neighborhood_width - 1U);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      if(!defects.defective(x, y)) {
        continue;
      }
      collect_same_phase_neighbors(measured, x, y, options.repair_radius,
                                   neighbors, &defects);
      if(neighbors.size() < options.min_neighbors) {
        astrocfa::CfaSample sample = cfa.sample_info(x, y);
        sample.valid = false;
        cfa.set_sample(x, y, sample);
        stats.unrepaired_pixels += 1;
        continue;
      }
      const float repaired = static_cast<float>(median_in_place(neighbors));
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = repaired,
                               .valid = true,
                               .clipped = false,
                           });
      stats.repaired_pixels += 1;
    }
  }
}

} // namespace

namespace astrocfa {

DefectMap::DefectMap(std::size_t width, std::size_t height)
    : width_(width), height_(height), flags_(width * height, 0U) {}

std::size_t DefectMap::offset(std::size_t x, std::size_t y) const {
  if(x >= width_ || y >= height_) {
    throw std::out_of_range("Defect map coordinate out of range");
  }
  return y * width_ + x;
}

SensorDefect DefectMap::at(std::size_t x, std::size_t y) const {
  return static_cast<SensorDefect>(flags_[offset(x, y)]);
}

bool DefectMap::has(std::size_t x, std::size_t y, SensorDefect defect) const {
  const auto flag = static_cast<std::uint8_t>(defect);
  return (flags_[offset(x, y)] & flag) != 0U;
}

bool DefectMap::defective(std::size_t x, std::size_t y) const {
  return flags_[offset(x, y)] != 0U;
}

void DefectMap::add(std::size_t x, std::size_t y, SensorDefect defect) {
  flags_[offset(x, y)] |= static_cast<std::uint8_t>(defect);
}

CalibrationResult calibrate_cfa(const CfaFrame &light, CalibrationInputs inputs) {
  require_compatible(light, inputs.bias, "Bias master");
  require_compatible(light, inputs.dark, "Dark master");
  require_compatible(light, inputs.flat, "Flat master");

  CalibrationResult result{
      .cfa = CfaFrame(light.width(), light.height(), light.pattern()),
      .defects = DefectMap(light.width(), light.height()),
  };
  result.defects = detect_sensor_defects(light, inputs.dark, inputs.flat,
                                         inputs.options.cosmetic, result.stats);
  compute_flat_phase_means(light, inputs.flat, result.stats.flat_phase_mean,
                           &result.defects);

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
    }
  }

  if(inputs.options.cosmetic.enabled) {
    repair_sensor_defects(result.cfa, result.defects, inputs.options.cosmetic,
                          result.stats);
  }

  const double usable = static_cast<double>(
      result.stats.samples - result.stats.invalid_samples - result.stats.clipped_samples);
  if(usable > 0.0) {
    result.stats.mean_before /= usable;
    result.stats.mean_bias_subtracted /= usable;
    result.stats.mean_dark_subtracted /= usable;
  }

  std::size_t usable_after = 0;
  for(std::size_t y = 0; y < result.cfa.height(); ++y) {
    for(std::size_t x = 0; x < result.cfa.width(); ++x) {
      const CfaSample sample = result.cfa.sample_info(x, y);
      if(sample.valid && !sample.clipped) {
        result.stats.mean_after += sample.value;
        usable_after += 1;
      }
    }
  }
  if(usable_after > 0) {
    result.stats.mean_after /= static_cast<double>(usable_after);
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
