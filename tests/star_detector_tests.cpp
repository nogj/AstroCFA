#include "astrocfa/star_detector.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void detects_bright_cfa_safe_star_candidate() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(8, 8, rggb);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{.value = 0.01F, .valid = true, .clipped = false});
    }
  }

  for(std::size_t y = 2; y < 4; ++y) {
    for(std::size_t x = 2; x < 4; ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = false});
    }
  }

  const astrocfa::StarDetectionStats stats = astrocfa::detect_star_candidates(frame, 2.0);
  require(stats.candidates == 1, "One bright synthetic star should be detected");
  require(stats.largest_area == 1, "The compact 2x2 star maps to one proxy sample");
  require(stats.brightest > 0.9F, "Brightest proxy sample should preserve stellar signal");
}

void ignores_clipped_blocks_when_estimating_background() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(4, 4, rggb);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{.value = 0.1F, .valid = true, .clipped = false});
    }
  }

  frame.set_sample(0, 0, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(1, 0, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(0, 1, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(1, 1, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});

  const astrocfa::StarDetectionStats stats = astrocfa::detect_star_candidates(frame, 3.0);
  require(stats.background_mean < 0.11, "Clipped block should not bias background upward");
}

} // namespace

int main() {
  try {
    detects_bright_cfa_safe_star_candidate();
    ignores_clipped_blocks_when_estimating_background();
  } catch(const std::exception &error) {
    std::cerr << "star_detector_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
