#include "astrocfa/synthetic_astro_scene.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace {

float clamp_scene(double value) {
  return static_cast<float>(std::clamp(value, 0.0, 1.25));
}

double sqr(double value) {
  return value * value;
}

double channel_value(astrocfa::RgbPixel pixel, int channel) {
  if(channel == 0) {
    return pixel.r;
  }
  if(channel == 2) {
    return pixel.b;
  }
  return pixel.g;
}

int channel_for_phase(astrocfa::CfaColor color) {
  switch(color) {
  case astrocfa::CfaColor::red:
    return 0;
  case astrocfa::CfaColor::green1:
  case astrocfa::CfaColor::green2:
    return 1;
  case astrocfa::CfaColor::blue:
    return 2;
  }
  return 1;
}

std::vector<astrocfa::SyntheticStar> make_stars(std::size_t width, std::size_t height) {
  return {
      {.x = width * 0.18, .y = height * 0.22, .flux = 0.72, .sigma = 0.48,
       .color = {.r = 1.0F, .g = 1.0F, .b = 1.0F}},
      {.x = width * 0.36 + 0.23, .y = height * 0.31 + 0.41, .flux = 0.58, .sigma = 0.62,
       .color = {.r = 0.88F, .g = 0.96F, .b = 1.0F}},
      {.x = width * 0.72 + 0.37, .y = height * 0.27 + 0.19, .flux = 0.85, .sigma = 0.52,
       .color = {.r = 1.0F, .g = 0.91F, .b = 0.76F}},
      {.x = width * 0.61 + 0.11, .y = height * 0.66 + 0.33, .flux = 1.20, .sigma = 0.42,
       .color = {.r = 1.0F, .g = 1.0F, .b = 1.0F}},
      {.x = width * 0.48 + 0.49, .y = height * 0.78 + 0.27, .flux = 0.34, .sigma = 0.95,
       .color = {.r = 0.72F, .g = 0.86F, .b = 1.0F}},
  };
}

} // namespace

namespace astrocfa {

SyntheticAstroScene make_synthetic_astro_scene(SyntheticAstroSceneOptions options) {
  RgbImage truth(options.width, options.height);
  std::vector<SyntheticStar> stars = make_stars(options.width, options.height);

  for(std::size_t y = 0; y < options.height; ++y) {
    const double fy = static_cast<double>(y) / static_cast<double>(options.height);
    for(std::size_t x = 0; x < options.width; ++x) {
      const double fx = static_cast<double>(x) / static_cast<double>(options.width);
      const double nebula =
          0.024 * std::exp(-(sqr((fx - 0.62) / 0.22) + sqr((fy - 0.55) / 0.16))) +
          0.010 * std::exp(-(sqr((fx - 0.34) / 0.18) + sqr((fy - 0.72) / 0.09)));
      const double background = 0.010 + 0.006 * fx + 0.004 * fy;
      RgbPixel pixel{
          .r = clamp_scene(background + 1.25 * nebula),
          .g = clamp_scene(background + 1.00 * nebula),
          .b = clamp_scene(background + 1.45 * nebula),
      };

      for(const SyntheticStar &star : stars) {
        const double dx = static_cast<double>(x) - star.x;
        const double dy = static_cast<double>(y) - star.y;
        const double psf = star.flux * std::exp(-0.5 * (dx * dx + dy * dy) /
                                                (star.sigma * star.sigma));
        pixel.r = clamp_scene(pixel.r + psf * star.color.r);
        pixel.g = clamp_scene(pixel.g + psf * star.color.g);
        pixel.b = clamp_scene(pixel.b + psf * star.color.b);
      }

      truth.set_pixel(x, y, pixel);
    }
  }

  CfaFrame cfa(options.width, options.height, options.pattern);
  std::mt19937 rng(options.seed);
  std::normal_distribution<double> gaussian(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> x_pick(0, options.width - 1U);
  std::uniform_int_distribution<std::size_t> y_pick(0, options.height - 1U);

  for(std::size_t y = 0; y < options.height; ++y) {
    for(std::size_t x = 0; x < options.width; ++x) {
      const int channel = channel_for_phase(options.pattern.at(x, y));
      double value = channel_value(truth.pixel(x, y), channel);
      if(options.add_noise) {
        const double sigma =
            std::sqrt(options.read_noise * options.read_noise +
                      std::max(0.0, value) * options.shot_noise_scale);
        value += gaussian(rng) * sigma;
      }
      CfaSample sample{
          .value = clamp_scene(value),
          .valid = true,
          .clipped = value >= 1.0,
      };
      cfa.set_sample(x, y, sample);
    }
  }

  if(options.add_hot_pixels) {
    const std::size_t hot_pixels = std::max<std::size_t>(1, options.width * options.height / 512U);
    for(std::size_t i = 0; i < hot_pixels; ++i) {
      const std::size_t x = x_pick(rng);
      const std::size_t y = y_pick(rng);
      cfa.set_sample(x, y, CfaSample{.value = 1.0F, .valid = true, .clipped = true});
    }
  }

  return SyntheticAstroScene{
      .truth = truth,
      .cfa = cfa,
      .stars = stars,
  };
}

} // namespace astrocfa
