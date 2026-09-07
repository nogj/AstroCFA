#include "astrocfa/star_detector.hpp"

#include "astrocfa/connected_components.hpp"
#include "astrocfa/luminance_proxy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace astrocfa {

StarDetectionStats detect_star_candidates(const CfaFrame &cfa, double sigma_threshold) {
  const LuminanceProxy luminance = build_cfa_luminance_proxy(cfa);
  const std::size_t proxy_width = luminance.width();
  const std::size_t proxy_height = luminance.height();
  StarDetectionStats stats;

  if(proxy_width == 0 || proxy_height == 0) {
    return stats;
  }

  std::vector<float> background_samples;
  background_samples.reserve(proxy_width * proxy_height);

  for(std::size_t by = 0; by < proxy_height; ++by) {
    for(std::size_t bx = 0; bx < proxy_width; ++bx) {
      const float value = luminance.value(bx, by);
      stats.brightest = std::max(stats.brightest, value);
      if(value > 0.0F) {
        background_samples.push_back(value);
      }
    }
  }

  if(background_samples.empty()) {
    return stats;
  }

  const RobustBackground background = robust_background(background_samples);
  stats.background_mean = background.location;
  stats.background_sigma = background.sigma;
  stats.threshold = stats.background_mean + sigma_threshold * stats.background_sigma;

  std::vector<std::uint8_t> mask(proxy_width * proxy_height, 0);
  for(std::size_t i = 0; i < mask.size(); ++i) {
    const std::size_t y = i / proxy_width;
    const std::size_t x = i - y * proxy_width;
    mask[i] = luminance.value(x, y) > stats.threshold ? 1 : 0;
  }

  const BinaryComponentStats components =
      analyze_binary_components(mask, proxy_width, proxy_height);
  stats.candidates = components.components;
  stats.largest_area = components.largest_area;
  return stats;
}

} // namespace astrocfa
