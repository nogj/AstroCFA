#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/frequency_cfa.hpp"

namespace astrocfa {

[[nodiscard]] RgbImage make_grayscale_map(std::size_t width, std::size_t height,
                                          float value);
[[nodiscard]] RgbImage make_frequency_alias_risk_map(const CfaFrame &cfa,
                                                     FrequencyCfaOptions options = {});
[[nodiscard]] RgbImage make_remosaic_residual_map(const CfaFrame &measured,
                                                  const RgbImage &reconstructed,
                                                  double display_scale = 250.0);

} // namespace astrocfa
