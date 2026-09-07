#pragma once

#include "astrocfa/cfa.hpp"

namespace astrocfa {

struct AstroPreviewOptions {
  double black_percentile = 0.10;
  double white_percentile = 99.80;
  double arcsinh_strength = 18.0;
};

[[nodiscard]] RgbImage make_astro_preview(const RgbImage &linear,
                                          AstroPreviewOptions options = {});

} // namespace astrocfa
