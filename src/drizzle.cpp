#include "astrocfa/drizzle.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astrocfa {

CfaDrizzleAccumulator::CfaDrizzleAccumulator(std::size_t input_width, std::size_t input_height,
                                             BayerPattern pattern, DrizzleOptions options)
    : input_width_(input_width), input_height_(input_height),
      output_width_(input_width * options.scale), output_height_(input_height * options.scale),
      pattern_(pattern), options_(options) {
  if(input_width == 0 || input_height == 0) {
    throw std::invalid_argument("Drizzle accumulator requires non-empty input dimensions");
  }
  if(options.scale == 0) {
    throw std::invalid_argument("Drizzle scale must be greater than zero");
  }
  if(options.drop_shrink <= 0.0 || options.drop_shrink > 1.0) {
    throw std::invalid_argument("Drizzle drop shrink must be in (0, 1]");
  }

  for(auto &plane : values_) {
    plane.assign(output_width_ * output_height_, 0.0);
  }
  for(auto &plane : weights_) {
    plane.assign(output_width_ * output_height_, 0.0);
  }
}

void CfaDrizzleAccumulator::add_frame(const CfaFrame &frame, SubpixelOffset offset,
                                      double frame_weight) {
  if(frame.width() != input_width_ || frame.height() != input_height_) {
    throw std::invalid_argument("Drizzle frame dimensions differ from accumulator");
  }
  if(frame_weight <= 0.0) {
    return;
  }

  const double scale = static_cast<double>(options_.scale);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      const CfaSample sample = frame.sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }

      const double target_x = (static_cast<double>(x) + offset.dx) * scale;
      const double target_y = (static_cast<double>(y) + offset.dy) * scale;
      deposit(frame.pattern().at(x, y), target_x, target_y, sample.value, frame_weight);
    }
  }
}

DrizzlePlaneSample CfaDrizzleAccumulator::plane_sample(CfaColor color, std::size_t x,
                                                       std::size_t y) const {
  const std::size_t sample_offset = offset(x, y);
  const std::size_t plane = plane_index(color);
  const double weight = weights_[plane][sample_offset];
  if(weight <= 0.0) {
    return {};
  }

  return DrizzlePlaneSample{
      .value = values_[plane][sample_offset] / weight,
      .weight = weight,
  };
}

DrizzleCoverageStats CfaDrizzleAccumulator::coverage(CfaColor color) const {
  const std::size_t plane = plane_index(color);
  DrizzleCoverageStats stats;
  stats.output_width = output_width_;
  stats.output_height = output_height_;
  stats.total_samples = output_width_ * output_height_;

  double weight_sum = 0.0;
  for(double weight : weights_[plane]) {
    if(weight > 0.0) {
      stats.covered_samples += 1;
      weight_sum += weight;
      stats.max_weight = std::max(stats.max_weight, weight);
    }
  }

  if(stats.covered_samples > 0) {
    stats.mean_weight = weight_sum / static_cast<double>(stats.covered_samples);
  }
  return stats;
}

std::size_t CfaDrizzleAccumulator::plane_index(CfaColor color) const {
  return static_cast<std::size_t>(color);
}

std::size_t CfaDrizzleAccumulator::offset(std::size_t x, std::size_t y) const {
  if(x >= output_width_ || y >= output_height_) {
    throw std::out_of_range("Drizzle output coordinate out of range");
  }
  return y * output_width_ + x;
}

void CfaDrizzleAccumulator::deposit(CfaColor color, double x, double y, float value,
                                    double weight) {
  const double radius = 0.5 * options_.drop_shrink * static_cast<double>(options_.scale);
  const double min_x = x - radius;
  const double max_x = x + radius;
  const double min_y = y - radius;
  const double max_y = y + radius;

  const int ix0 = static_cast<int>(std::floor(min_x));
  const int ix1 = static_cast<int>(std::ceil(max_x));
  const int iy0 = static_cast<int>(std::floor(min_y));
  const int iy1 = static_cast<int>(std::ceil(max_y));
  const std::size_t plane = plane_index(color);

  double normalization = 0.0;
  for(int iy = iy0; iy <= iy1; ++iy) {
    for(int ix = ix0; ix <= ix1; ++ix) {
      if(ix < 0 || iy < 0 || ix >= static_cast<int>(output_width_) ||
         iy >= static_cast<int>(output_height_)) {
        continue;
      }

      const double center_x = static_cast<double>(ix) + 0.5;
      const double center_y = static_cast<double>(iy) + 0.5;
      const double wx = std::max(0.0, radius - std::abs(center_x - x));
      const double wy = std::max(0.0, radius - std::abs(center_y - y));
      normalization += wx * wy;
    }
  }

  if(normalization <= 0.0) {
    return;
  }

  for(int iy = iy0; iy <= iy1; ++iy) {
    for(int ix = ix0; ix <= ix1; ++ix) {
      if(ix < 0 || iy < 0 || ix >= static_cast<int>(output_width_) ||
         iy >= static_cast<int>(output_height_)) {
        continue;
      }

      const double center_x = static_cast<double>(ix) + 0.5;
      const double center_y = static_cast<double>(iy) + 0.5;
      const double wx = std::max(0.0, radius - std::abs(center_x - x));
      const double wy = std::max(0.0, radius - std::abs(center_y - y));
      const double sample_weight = (wx * wy / normalization) * weight;
      if(sample_weight <= 0.0) {
        continue;
      }
      const std::size_t target = offset(static_cast<std::size_t>(ix), static_cast<std::size_t>(iy));
      values_[plane][target] += static_cast<double>(value) * sample_weight;
      weights_[plane][target] += sample_weight;
    }
  }
}

} // namespace astrocfa
