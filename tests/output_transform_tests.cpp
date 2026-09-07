#include "astrocfa/output_transform.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

bool near(double left, double right, double tolerance = 1.0e-5) {
  return std::abs(left - right) <= tolerance;
}

astrocfa::RgbImage grayscale_ramp() {
  astrocfa::RgbImage image(101, 1);
  for(std::size_t x = 0; x < image.width(); ++x) {
    const float value = static_cast<float>(x) / 100.0F;
    image.set_pixel(x, 0, {.r = value, .g = value, .b = value});
  }
  return image;
}

void legacy_preview_is_monotonic() {
  const astrocfa::RgbImage preview =
      astrocfa::make_astro_preview(grayscale_ramp());
  require(preview.width() == 101 && preview.height() == 1,
          "Preview should preserve dimensions");
  for(std::size_t x = 1; x < preview.width(); ++x) {
    require(preview.pixel(x - 1U, 0).g <= preview.pixel(x, 0).g,
            "Astro preview should be monotonic");
  }
}

void ghs_is_normalized_and_monotonic() {
  const astrocfa::AstroToneResult result = astrocfa::apply_astro_tone(
      grayscale_ramp(), astrocfa::AstroToneOptions{
                            .curve =
                                astrocfa::ToneCurve::generalized_hyperbolic,
                            .auto_levels = false,
                            .stretch_factor = 3.0,
                            .local_intensity = 8.0,
                            .symmetry_point = 0.08,
                            .shadow_protection = 0.02,
                            .highlight_protection = 0.80,
                        });
  require(near(result.image.pixel(0, 0).g, 0.0),
          "GHS should map zero to zero");
  require(near(result.image.pixel(100, 0).g, 1.0),
          "GHS should map one to one");
  for(std::size_t x = 1; x < result.image.width(); ++x) {
    require(result.image.pixel(x - 1U, 0).g <= result.image.pixel(x, 0).g,
            "GHS should remain monotonic");
  }
}

void luminance_mode_preserves_color_ratios() {
  astrocfa::RgbImage image(1, 1);
  image.set_pixel(0, 0, {.r = 0.10F, .g = 0.20F, .b = 0.30F});
  const astrocfa::RgbPixel output =
      astrocfa::apply_astro_tone(
          image, astrocfa::AstroToneOptions{
                     .curve = astrocfa::ToneCurve::arcsinh,
                     .auto_levels = false,
                     .stretch_factor = 2.0,
                 })
          .image.pixel(0, 0);
  require(near(output.r / output.g, 0.5) &&
              near(output.b / output.g, 1.5),
          "Luminance stretch should preserve RGB ratios when gamut allows");
}

void linear_curve_applies_exposure_without_gamma() {
  astrocfa::RgbImage image(1, 1);
  image.set_pixel(0, 0, {.r = 0.25F, .g = 0.25F, .b = 0.25F});
  const astrocfa::AstroToneResult result = astrocfa::apply_astro_tone(
      image, astrocfa::AstroToneOptions{
                 .curve = astrocfa::ToneCurve::linear,
                 .exposure_ev = 1.0,
                 .auto_levels = false,
             });
  require(near(result.image.pixel(0, 0).g, 0.5),
          "Linear tone should apply exposure in scene-linear space");
}

void highlight_segment_reserves_bright_star_contrast() {
  const astrocfa::RgbImage ramp = grayscale_ramp();
  const astrocfa::AstroToneOptions base{
      .curve = astrocfa::ToneCurve::generalized_hyperbolic,
      .auto_levels = false,
      .stretch_factor = 3.0,
      .local_intensity = 8.0,
      .symmetry_point = 0.08,
  };
  astrocfa::AstroToneOptions protected_options = base;
  protected_options.highlight_protection = 0.80;
  const astrocfa::RgbImage unprotected =
      astrocfa::apply_astro_tone(ramp, base).image;
  const astrocfa::RgbImage protected_image =
      astrocfa::apply_astro_tone(ramp, protected_options).image;
  const double unprotected_contrast =
      unprotected.pixel(95, 0).g - unprotected.pixel(90, 0).g;
  const double protected_contrast =
      protected_image.pixel(95, 0).g - protected_image.pixel(90, 0).g;
  require(protected_contrast > unprotected_contrast,
          "Highlight protection should reserve more contrast for bright stars");
}

} // namespace

int main() {
  try {
    legacy_preview_is_monotonic();
    ghs_is_normalized_and_monotonic();
    luminance_mode_preserves_color_ratios();
    linear_curve_applies_exposure_without_gamma();
    highlight_segment_reserves_bright_star_contrast();
  } catch(const std::exception &error) {
    std::cerr << "output_transform_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
