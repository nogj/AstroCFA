#include "astrocfa/registration.hpp"

#include "astrocfa/luminance_proxy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

struct ThresholdedProxy {
  astrocfa::LuminanceProxy proxy;
  double threshold = 0.0;
};

ThresholdedProxy make_thresholded_proxy(const astrocfa::CfaFrame &cfa,
                                        double sigma_threshold) {
  astrocfa::LuminanceProxy proxy = astrocfa::build_cfa_luminance_proxy(cfa);
  std::vector<float> samples;
  samples.reserve(proxy.width() * proxy.height());

  for(std::size_t y = 0; y < proxy.height(); ++y) {
    for(std::size_t x = 0; x < proxy.width(); ++x) {
      const float value = proxy.value(x, y);
      if(value > 0.0F) {
        samples.push_back(value);
      }
    }
  }

  const astrocfa::RobustBackground background = astrocfa::robust_background(samples);
  return ThresholdedProxy{
      .proxy = proxy,
      .threshold = background.location + sigma_threshold * background.sigma,
  };
}

double signal_above_threshold(const astrocfa::LuminanceProxy &proxy, std::size_t x,
                              std::size_t y, double threshold) {
  return std::max(0.0, static_cast<double>(proxy.value(x, y)) - threshold);
}

} // namespace

namespace astrocfa {

RegistrationResult estimate_integer_star_translation(const CfaFrame &reference,
                                                     const CfaFrame &moving,
                                                     RegistrationOptions options) {
  if(reference.width() != moving.width() || reference.height() != moving.height()) {
    throw std::invalid_argument("Registration frames must have matching dimensions");
  }

  const ThresholdedProxy ref = make_thresholded_proxy(reference, options.sigma_threshold);
  const ThresholdedProxy mov = make_thresholded_proxy(moving, options.sigma_threshold);
  if(ref.proxy.width() != mov.proxy.width() || ref.proxy.height() != mov.proxy.height()) {
    throw std::invalid_argument("Registration proxies must have matching dimensions");
  }

  const int max_shift_proxy = std::max(0, options.max_shift_pixels / 2);
  const std::size_t stride = std::max<std::size_t>(1, options.sample_stride);
  RegistrationResult best;
  best.score = -std::numeric_limits<double>::infinity();

  const auto width = static_cast<int>(ref.proxy.width());
  const auto height = static_cast<int>(ref.proxy.height());

  for(int dy = -max_shift_proxy; dy <= max_shift_proxy; ++dy) {
    for(int dx = -max_shift_proxy; dx <= max_shift_proxy; ++dx) {
      double score = 0.0;
      std::size_t matched = 0;

      for(int y = 0; y < height; y += static_cast<int>(stride)) {
        const int my = y + dy;
        if(my < 0 || my >= height) {
          continue;
        }

        for(int x = 0; x < width; x += static_cast<int>(stride)) {
          const int mx = x + dx;
          if(mx < 0 || mx >= width) {
            continue;
          }

          const double ref_signal = signal_above_threshold(
              ref.proxy, static_cast<std::size_t>(x), static_cast<std::size_t>(y),
              ref.threshold);
          const double mov_signal = signal_above_threshold(
              mov.proxy, static_cast<std::size_t>(mx), static_cast<std::size_t>(my),
              mov.threshold);
          if(ref_signal <= 0.0 && mov_signal <= 0.0) {
            continue;
          }

          score += ref_signal * mov_signal;
          matched += 1;
        }
      }

      if(score > best.score) {
        best.score = score;
        best.proxy_dx = -dx;
        best.proxy_dy = -dy;
        best.offset = SubpixelOffset{
            .dx = static_cast<double>(best.proxy_dx * 2),
            .dy = static_cast<double>(best.proxy_dy * 2),
        };
        best.matched_samples = matched;
      }
    }
  }

  if(!std::isfinite(best.score) || best.score <= 0.0) {
    best.score = 0.0;
    best.proxy_dx = 0;
    best.proxy_dy = 0;
    best.offset = SubpixelOffset{};
    best.matched_samples = 0;
  }
  return best;
}

} // namespace astrocfa
