#pragma once

#include "astrocfa/cfa.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace astrocfa {

struct SubpixelOffset {
  double dx = 0.0;
  double dy = 0.0;
};

struct DrizzleOptions {
  std::size_t scale = 2;
  double drop_shrink = 0.7;
};

struct DrizzlePlaneSample {
  double value = 0.0;
  double weight = 0.0;
};

struct DrizzleCoverageStats {
  std::size_t output_width = 0;
  std::size_t output_height = 0;
  std::size_t covered_samples = 0;
  std::size_t total_samples = 0;
  double mean_weight = 0.0;
  double max_weight = 0.0;
};

class CfaDrizzleAccumulator {
public:
  CfaDrizzleAccumulator(std::size_t input_width, std::size_t input_height,
                        BayerPattern pattern, DrizzleOptions options);

  void add_frame(const CfaFrame &frame, SubpixelOffset offset, double frame_weight = 1.0);

  [[nodiscard]] std::size_t output_width() const { return output_width_; }
  [[nodiscard]] std::size_t output_height() const { return output_height_; }
  [[nodiscard]] const BayerPattern &pattern() const { return pattern_; }
  [[nodiscard]] DrizzleOptions options() const { return options_; }

  [[nodiscard]] DrizzlePlaneSample plane_sample(CfaColor color, std::size_t x,
                                                std::size_t y) const;
  [[nodiscard]] DrizzleCoverageStats coverage(CfaColor color) const;

private:
  [[nodiscard]] std::size_t plane_index(CfaColor color) const;
  [[nodiscard]] std::size_t offset(std::size_t x, std::size_t y) const;
  void deposit(CfaColor color, double x, double y, float value, double weight);

  std::size_t input_width_ = 0;
  std::size_t input_height_ = 0;
  std::size_t output_width_ = 0;
  std::size_t output_height_ = 0;
  BayerPattern pattern_;
  DrizzleOptions options_;
  std::array<std::vector<double>, 4> values_;
  std::array<std::vector<double>, 4> weights_;
};

} // namespace astrocfa

