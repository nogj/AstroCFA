#include "astrocfa/noise_model.hpp"

#include <algorithm>
#include <cmath>

namespace astrocfa {

NoiseEstimate estimate_noise(double normalized_signal, const NoiseModel &model) {
  const double signal = std::max(0.0, normalized_signal);
  const double variance =
      std::max(model.floor, model.read_noise * model.read_noise +
                                model.shot_noise_scale * signal +
                                model.quantization_noise * model.quantization_noise);

  NoiseEstimate estimate;
  estimate.variance = variance;
  estimate.sigma = std::sqrt(variance);
  estimate.weight = 1.0 / variance;
  return estimate;
}

} // namespace astrocfa

