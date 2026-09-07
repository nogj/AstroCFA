#include "astrocfa/diagnostic_maps.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

astrocfa::RgbPixel gray(float value) {
  value = std::clamp(value, 0.0F, 1.0F);
  return astrocfa::RgbPixel{.r = value, .g = value, .b = value};
}

} // namespace

namespace astrocfa {

RgbImage make_grayscale_map(std::size_t width, std::size_t height, float value) {
  RgbImage image(width, height);
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      image.set_pixel(x, y, gray(value));
    }
  }
  return image;
}

RgbImage make_frequency_alias_risk_map(const CfaFrame &cfa,
                                       FrequencyCfaOptions options) {
  const FrequencyRiskMap risk = build_frequency_risk_map(cfa, options);
  RgbImage image(cfa.width(), cfa.height());
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      image.set_pixel(x, y, gray(static_cast<float>(risk.pixel_risk(x, y))));
    }
  }
  return image;
}

RgbImage make_remosaic_residual_map(const CfaFrame &measured,
                                    const RgbImage &reconstructed,
                                    double display_scale) {
  if(measured.width() != reconstructed.width() ||
     measured.height() != reconstructed.height()) {
    throw std::invalid_argument("Measured CFA and reconstructed RGB dimensions differ");
  }
  if(display_scale <= 0.0) {
    throw std::invalid_argument("Residual display scale must be positive");
  }

  RgbImage image(measured.width(), measured.height());
  for(std::size_t y = 0; y < measured.height(); ++y) {
    for(std::size_t x = 0; x < measured.width(); ++x) {
      const CfaSample sample = measured.sample_info(x, y);
      if(!sample.valid) {
        image.set_pixel(x, y, RgbPixel{.r = 0.0F, .g = 0.0F, .b = 1.0F});
        continue;
      }
      if(sample.clipped) {
        image.set_pixel(x, y, RgbPixel{.r = 1.0F, .g = 0.0F, .b = 0.0F});
        continue;
      }

      const float predicted =
          remosaic_pixel(reconstructed.pixel(x, y), measured.pattern().at(x, y));
      const float residual =
          static_cast<float>(std::abs(static_cast<double>(predicted) - sample.value) *
                             display_scale);
      image.set_pixel(x, y, gray(residual));
    }
  }
  return image;
}

} // namespace astrocfa
