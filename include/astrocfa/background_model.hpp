#pragma once

#include "astrocfa/cfa.hpp"

#include <array>
#include <cstddef>

namespace astrocfa {

struct BackgroundModelOptions {
  std::size_t tile_size = 64;
  std::size_t polynomial_degree = 2;
  double sample_quantile = 0.20;
  std::size_t robust_iterations = 6;
  double positive_rejection_sigma = 2.5;
  double negative_rejection_sigma = 4.0;
  bool neutralize = false;
};

struct BackgroundModelStats {
  std::size_t effective_tile_size = 0;
  std::size_t tile_samples = 0;
  std::size_t downweighted_samples = 0;
  std::array<double, 3> preserved_level = {0.0, 0.0, 0.0};
  std::array<double, 3> gradient_peak_to_peak = {0.0, 0.0, 0.0};
  double robust_residual_sigma = 0.0;
};

struct BackgroundModelResult {
  RgbImage corrected;
  RgbImage background;
  BackgroundModelStats stats;
};

[[nodiscard]] BackgroundModelResult model_astro_background(
    const RgbImage &linear_camera_rgb, BackgroundModelOptions options = {});

} // namespace astrocfa
