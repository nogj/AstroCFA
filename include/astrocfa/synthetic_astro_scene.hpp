#pragma once

#include "astrocfa/cfa.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>

namespace astrocfa {

struct SyntheticStar {
  double x = 0.0;
  double y = 0.0;
  double flux = 0.0;
  double sigma = 1.0;
  RgbPixel color{.r = 1.0F, .g = 1.0F, .b = 1.0F};
};

struct SyntheticAstroSceneOptions {
  std::size_t width = 128;
  std::size_t height = 96;
  std::uint32_t seed = 42;
  bool add_noise = true;
  bool add_hot_pixels = true;
  double read_noise = 0.0025;
  double shot_noise_scale = 0.0018;
  BayerPattern pattern;
};

struct SyntheticAstroScene {
  RgbImage truth;
  CfaFrame cfa;
  std::vector<SyntheticStar> stars;
};

[[nodiscard]] SyntheticAstroScene
make_synthetic_astro_scene(SyntheticAstroSceneOptions options = {});

} // namespace astrocfa
