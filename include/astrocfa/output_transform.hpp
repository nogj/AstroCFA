#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>

namespace astrocfa {

enum class ToneCurve {
  linear,
  arcsinh,
  generalized_hyperbolic,
};

struct AstroToneOptions {
  ToneCurve curve = ToneCurve::arcsinh;
  double exposure_ev = 0.0;
  bool auto_levels = true;
  double black_percentile = 0.10;
  double white_percentile = 99.80;
  double black_point = 0.0;
  double white_point = 1.0;
  double stretch_factor = 2.944438979;
  double local_intensity = 0.0;
  double symmetry_point = 0.0;
  double shadow_protection = 0.0;
  double highlight_protection = 1.0;
  double saturation = 1.0;
};

struct AstroToneStats {
  double black_point = 0.0;
  double white_point = 1.0;
  std::size_t shadow_clipped_pixels = 0;
  std::size_t highlight_clipped_pixels = 0;
  std::size_t gamut_compressed_pixels = 0;
  double maximum_input_luminance = 0.0;
  double maximum_output_component = 0.0;
};

struct AstroToneResult {
  RgbImage image;
  AstroToneStats stats;
};

struct AstroPreviewOptions {
  double black_percentile = 0.10;
  double white_percentile = 99.80;
  double arcsinh_strength = 18.0;
};

[[nodiscard]] RgbImage make_astro_preview(const RgbImage &linear,
                                          AstroPreviewOptions options = {});
[[nodiscard]] AstroToneResult apply_astro_tone(
    const RgbImage &linear, AstroToneOptions options = {});

} // namespace astrocfa
