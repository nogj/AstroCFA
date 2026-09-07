#include "astrocfa/luminance_proxy.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astrocfa {

LuminanceProxy::LuminanceProxy(std::size_t width, std::size_t height)
    : width_(width), height_(height), values_(width * height, 0.0F) {}

float LuminanceProxy::value(std::size_t x, std::size_t y) const {
  return values_[offset(x, y)];
}

void LuminanceProxy::set_value(std::size_t x, std::size_t y, float value) {
  values_[offset(x, y)] = value;
}

std::size_t LuminanceProxy::offset(std::size_t x, std::size_t y) const {
  if(x >= width_ || y >= height_) {
    throw std::out_of_range("Luminance proxy coordinate out of range");
  }
  return y * width_ + x;
}

double median(std::vector<float> values) {
  if(values.empty()) {
    return 0.0;
  }

  const auto middle = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2U);
  std::nth_element(values.begin(), middle, values.end());
  double result = *middle;

  if(values.size() % 2U == 0) {
    const auto lower = std::max_element(values.begin(), middle);
    result = (result + *lower) * 0.5;
  }

  return result;
}

RobustBackground robust_background(const std::vector<float> &samples) {
  RobustBackground background;
  background.location = median(samples);

  std::vector<float> absolute_deviations;
  absolute_deviations.reserve(samples.size());
  for(float sample : samples) {
    absolute_deviations.push_back(static_cast<float>(std::abs(sample - background.location)));
  }

  background.sigma = 1.4826 * median(absolute_deviations);
  return background;
}

LuminanceProxy build_cfa_luminance_proxy(const CfaFrame &cfa) {
  LuminanceProxy proxy(cfa.width() / 2U, cfa.height() / 2U);

  for(std::size_t by = 0; by < proxy.height(); ++by) {
    for(std::size_t bx = 0; bx < proxy.width(); ++bx) {
      double block_sum = 0.0;
      std::size_t block_samples = 0;

      for(std::size_t dy = 0; dy < 2; ++dy) {
        for(std::size_t dx = 0; dx < 2; ++dx) {
          const CfaSample sample = cfa.sample_info(bx * 2U + dx, by * 2U + dy);
          if(!sample.valid || sample.clipped) {
            continue;
          }
          block_sum += sample.value;
          block_samples += 1;
        }
      }

      proxy.set_value(
          bx, by,
          block_samples > 0 ? static_cast<float>(block_sum / static_cast<double>(block_samples))
                            : 0.0F);
    }
  }

  return proxy;
}

} // namespace astrocfa

