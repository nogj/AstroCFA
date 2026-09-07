#pragma once

#include <cstddef>

namespace astrocfa {

struct NoiseModel {
  double read_noise = 0.003;
  double shot_noise_scale = 0.001;
  double quantization_noise = 0.0;
  double floor = 1.0e-6;
};

struct NoiseEstimate {
  double variance = 0.0;
  double sigma = 0.0;
  double weight = 0.0;
};

struct NoiseWeightedResidual {
  double chi_square = 0.0;
  double reduced_chi_square = 0.0;
  double mean_normalized_absolute = 0.0;
  double max_normalized_absolute = 0.0;
  std::size_t samples = 0;
};

[[nodiscard]] NoiseEstimate estimate_noise(double normalized_signal,
                                           const NoiseModel &model);

} // namespace astrocfa

