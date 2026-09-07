#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/synthetic_astro_scene.hpp"

#include <cstddef>

namespace astrocfa {

struct ReconstructionMetrics {
  double rgb_mae = 0.0;
  double rgb_rmse = 0.0;
  double max_abs_error = 0.0;
  double chroma_mae = 0.0;
  double star_false_color = 0.0;
  double star_luma_rmse = 0.0;
  double star_flux_relative_error = 0.0;
  double star_flux_relative_bias = 0.0;
  double star_fwhm_relative_error = 0.0;
  double star_elongation_error = 0.0;
  double cfa_residual_mae = 0.0;
  std::size_t samples = 0;
  std::size_t star_samples = 0;
  std::size_t measured_stars = 0;
};

[[nodiscard]] ReconstructionMetrics
measure_reconstruction(const RgbImage &truth, const RgbImage &reconstructed,
                       const CfaFrame &measured,
                       const std::vector<SyntheticStar> &stars);

} // namespace astrocfa
