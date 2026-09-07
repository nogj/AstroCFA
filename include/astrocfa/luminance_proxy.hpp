#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>
#include <vector>

namespace astrocfa {

class LuminanceProxy {
public:
  LuminanceProxy(std::size_t width, std::size_t height);

  [[nodiscard]] std::size_t width() const { return width_; }
  [[nodiscard]] std::size_t height() const { return height_; }
  [[nodiscard]] float value(std::size_t x, std::size_t y) const;
  void set_value(std::size_t x, std::size_t y, float value);

private:
  [[nodiscard]] std::size_t offset(std::size_t x, std::size_t y) const;

  std::size_t width_ = 0;
  std::size_t height_ = 0;
  std::vector<float> values_;
};

struct RobustBackground {
  double location = 0.0;
  double sigma = 0.0;
};

[[nodiscard]] double median(std::vector<float> values);
[[nodiscard]] RobustBackground robust_background(const std::vector<float> &samples);
[[nodiscard]] LuminanceProxy build_cfa_luminance_proxy(const CfaFrame &cfa);

} // namespace astrocfa

