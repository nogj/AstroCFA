#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "astrocfa/noise_model.hpp"

namespace astrocfa {

enum class CfaColor : int {
  red = 0,
  green1 = 1,
  blue = 2,
  green2 = 3,
};

struct BayerPattern {
  std::array<CfaColor, 4> phase = {
      CfaColor::red,
      CfaColor::green1,
      CfaColor::green2,
      CfaColor::blue,
  };

  [[nodiscard]] CfaColor at(std::size_t x, std::size_t y) const {
    return phase[(y % 2U) * 2U + (x % 2U)];
  }
};

struct RgbPixel {
  float r = 0.0F;
  float g = 0.0F;
  float b = 0.0F;
};

struct CfaSample {
  float value = 0.0F;
  bool valid = true;
  bool clipped = false;
};

class CfaFrame {
public:
  CfaFrame(std::size_t width, std::size_t height, BayerPattern pattern);

  [[nodiscard]] std::size_t width() const { return width_; }
  [[nodiscard]] std::size_t height() const { return height_; }
  [[nodiscard]] const BayerPattern &pattern() const { return pattern_; }

  [[nodiscard]] float sample(std::size_t x, std::size_t y) const;
  [[nodiscard]] CfaSample sample_info(std::size_t x, std::size_t y) const;
  void set_sample(std::size_t x, std::size_t y, float value);
  void set_sample(std::size_t x, std::size_t y, CfaSample sample);
  void set_valid(std::size_t x, std::size_t y, bool valid);
  void set_clipped(std::size_t x, std::size_t y, bool clipped);

private:
  [[nodiscard]] std::size_t offset(std::size_t x, std::size_t y) const;

  std::size_t width_ = 0;
  std::size_t height_ = 0;
  BayerPattern pattern_;
  std::vector<CfaSample> samples_;
};

class RgbImage {
public:
  RgbImage(std::size_t width, std::size_t height);

  [[nodiscard]] std::size_t width() const { return width_; }
  [[nodiscard]] std::size_t height() const { return height_; }

  [[nodiscard]] RgbPixel pixel(std::size_t x, std::size_t y) const;
  void set_pixel(std::size_t x, std::size_t y, RgbPixel value);

private:
  [[nodiscard]] std::size_t offset(std::size_t x, std::size_t y) const;

  std::size_t width_ = 0;
  std::size_t height_ = 0;
  std::vector<RgbPixel> pixels_;
};

struct RemosaicResidual {
  double mean_absolute = 0.0;
  double root_mean_square = 0.0;
  float maximum_absolute = 0.0F;
  std::size_t samples = 0;
};

[[nodiscard]] float remosaic_pixel(RgbPixel pixel, CfaColor color);
[[nodiscard]] CfaFrame remosaic(const RgbImage &image, BayerPattern pattern);
[[nodiscard]] RemosaicResidual compute_remosaic_residual(const CfaFrame &measured,
                                                         const RgbImage &reconstructed);
[[nodiscard]] NoiseWeightedResidual
compute_noise_weighted_remosaic_residual(const CfaFrame &measured,
                                         const RgbImage &reconstructed,
                                         const NoiseModel &noise_model);

} // namespace astrocfa
