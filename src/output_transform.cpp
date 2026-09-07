#include "astrocfa/output_transform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

float luminance(astrocfa::RgbPixel pixel) {
  return 0.2126F * pixel.r + 0.7152F * pixel.g + 0.0722F * pixel.b;
}

double percentile(std::vector<float> values, double p) {
  if(values.empty()) {
    return 0.0;
  }

  p = std::clamp(p, 0.0, 100.0);
  const double index = (p / 100.0) * static_cast<double>(values.size() - 1U);
  const auto lower = static_cast<std::size_t>(std::floor(index));
  const auto upper = static_cast<std::size_t>(std::ceil(index));
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(lower),
                   values.end());
  const double lower_value = values[lower];
  if(upper == lower) {
    return lower_value;
  }
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(upper),
                   values.end());
  const double upper_value = values[upper];
  const double fraction = index - static_cast<double>(lower);
  return lower_value * (1.0 - fraction) + upper_value * fraction;
}

float map_channel(float value, double black, double white, double strength,
                  double normalizer) {
  if(white <= black) {
    return 0.0F;
  }

  const double normalized =
      std::clamp((static_cast<double>(value) - black) / (white - black), 0.0, 1.0);
  return static_cast<float>(std::asinh(strength * normalized) / normalizer);
}

} // namespace

namespace astrocfa {

RgbImage make_astro_preview(const RgbImage &linear, AstroPreviewOptions options) {
  if(options.arcsinh_strength <= 0.0) {
    throw std::invalid_argument("Astro preview arcsinh strength must be positive");
  }

  std::vector<float> luminance_samples;
  luminance_samples.reserve(linear.width() * linear.height());
  for(std::size_t y = 0; y < linear.height(); ++y) {
    for(std::size_t x = 0; x < linear.width(); ++x) {
      const float value = luminance(linear.pixel(x, y));
      if(std::isfinite(value)) {
        luminance_samples.push_back(value);
      }
    }
  }

  const double black = percentile(luminance_samples, options.black_percentile);
  const double white = percentile(luminance_samples, options.white_percentile);
  const double normalizer = std::asinh(options.arcsinh_strength);

  RgbImage preview(linear.width(), linear.height());
  for(std::size_t y = 0; y < linear.height(); ++y) {
    for(std::size_t x = 0; x < linear.width(); ++x) {
      const RgbPixel pixel = linear.pixel(x, y);
      preview.set_pixel(x, y,
                        RgbPixel{
                            .r = map_channel(pixel.r, black, white,
                                             options.arcsinh_strength, normalizer),
                            .g = map_channel(pixel.g, black, white,
                                             options.arcsinh_strength, normalizer),
                            .b = map_channel(pixel.b, black, white,
                                             options.arcsinh_strength, normalizer),
                        });
    }
  }

  return preview;
}

} // namespace astrocfa
