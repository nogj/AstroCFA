#include "astrocfa/background_model.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

astrocfa::RgbPixel synthetic_background(std::size_t x, std::size_t y,
                                        std::size_t width,
                                        std::size_t height) {
  const double nx = 2.0 * x / static_cast<double>(width - 1U) - 1.0;
  const double ny = 2.0 * y / static_cast<double>(height - 1U) - 1.0;
  return {
      .r = static_cast<float>(0.045 + 0.018 * nx + 0.009 * ny + 0.004 * nx * ny),
      .g = static_cast<float>(0.040 - 0.012 * nx + 0.007 * ny),
      .b = static_cast<float>(0.052 + 0.006 * nx - 0.014 * ny),
  };
}

astrocfa::RgbImage synthetic_field(bool include_sources) {
  constexpr std::size_t width = 128;
  constexpr std::size_t height = 96;
  astrocfa::RgbImage image(width, height);
  const double stars[][3] = {
      {20.0, 18.0, 0.40}, {42.0, 70.0, 0.55}, {92.0, 22.0, 0.35},
      {108.0, 72.0, 0.60}, {65.0, 46.0, 0.45},
  };
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      astrocfa::RgbPixel pixel = synthetic_background(x, y, width, height);
      if(include_sources) {
        double source = 0.0;
        for(const auto &star : stars) {
          const double dx = static_cast<double>(x) - star[0];
          const double dy = static_cast<double>(y) - star[1];
          source += star[2] * std::exp(-0.5 * (dx * dx + dy * dy) / 2.25);
        }
        const double nebula_dx = static_cast<double>(x) - 68.0;
        const double nebula_dy = static_cast<double>(y) - 49.0;
        const double nebula =
            0.075 * std::exp(-0.5 * (nebula_dx * nebula_dx / 500.0 +
                                     nebula_dy * nebula_dy / 220.0));
        pixel.r += static_cast<float>(source + 1.15 * nebula);
        pixel.g += static_cast<float>(source + 0.90 * nebula);
        pixel.b += static_cast<float>(source + 1.30 * nebula);
      }
      image.set_pixel(x, y, pixel);
    }
  }
  return image;
}

double channel_range_at_corners(const astrocfa::RgbImage &image,
                                std::size_t channel) {
  const std::size_t xs[] = {4U, image.width() - 5U};
  const std::size_t ys[] = {4U, image.height() - 5U};
  double minimum = 1.0e9;
  double maximum = -1.0e9;
  for(std::size_t y : ys) {
    for(std::size_t x : xs) {
      const astrocfa::RgbPixel pixel = image.pixel(x, y);
      const double value = channel == 0 ? pixel.r : channel == 1 ? pixel.g : pixel.b;
      minimum = std::min(minimum, value);
      maximum = std::max(maximum, value);
    }
  }
  return maximum - minimum;
}

void removes_gradient_without_modeling_extended_signal() {
  const astrocfa::RgbImage field = synthetic_field(true);
  const astrocfa::BackgroundModelResult result =
      astrocfa::model_astro_background(
          field, astrocfa::BackgroundModelOptions{.tile_size = 16});
  for(std::size_t channel = 0; channel < 3; ++channel) {
    require(channel_range_at_corners(result.corrected, channel) < 0.008,
            "Corrected corners should have a nearly uniform background");
  }

  const astrocfa::RgbPixel modeled = result.background.pixel(68, 49);
  const astrocfa::RgbPixel truth = synthetic_background(68, 49, 128, 96);
  require(std::abs(modeled.r - truth.r) < 0.018 &&
              std::abs(modeled.g - truth.g) < 0.018 &&
              std::abs(modeled.b - truth.b) < 0.018,
          "Robust lower-envelope fit should not absorb the central nebula");
  require(result.stats.downweighted_samples > 0,
          "Source-contaminated tiles should be downweighted");
}

void preserves_additive_star_signal() {
  const astrocfa::BackgroundModelOptions options{.tile_size = 16};
  const astrocfa::BackgroundModelResult sourced =
      astrocfa::model_astro_background(synthetic_field(true), options);
  const astrocfa::BackgroundModelResult empty =
      astrocfa::model_astro_background(synthetic_field(false), options);
  const double recovered_star = sourced.corrected.pixel(20, 18).g -
                                empty.corrected.pixel(20, 18).g;
  require(std::abs(recovered_star - 0.40) < 0.025,
          "Background correction should preserve compact stellar signal");
}

void neutralization_uses_one_preserved_sky_level() {
  const astrocfa::BackgroundModelResult result =
      astrocfa::model_astro_background(
          synthetic_field(false),
          astrocfa::BackgroundModelOptions{.tile_size = 16, .neutralize = true});
  require(result.stats.preserved_level[0] == result.stats.preserved_level[1] &&
              result.stats.preserved_level[1] == result.stats.preserved_level[2],
          "Neutral mode should use a common sky level across RGB");
}

void adapts_tile_size_for_small_images() {
  const astrocfa::BackgroundModelResult result =
      astrocfa::model_astro_background(synthetic_field(false));
  require(result.stats.effective_tile_size < 64,
          "Small images should automatically refine the background tile grid");
  require(result.stats.tile_samples >= 8,
          "Adaptive grid should retain enough samples for a quadratic fit");
}

} // namespace

int main() {
  try {
    removes_gradient_without_modeling_extended_signal();
    preserves_additive_star_signal();
    neutralization_uses_one_preserved_sky_level();
    adapts_tile_size_for_small_images();
  } catch(const std::exception &error) {
    std::cerr << "background_model_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
