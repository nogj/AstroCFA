#include "astrocfa/registration.hpp"

#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

astrocfa::CfaFrame synthetic_star_field(std::size_t width, std::size_t height,
                                        int star_x, int star_y) {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(width, height, rggb);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(x, y,
                       astrocfa::CfaSample{.value = 0.01F, .valid = true, .clipped = false});
    }
  }

  for(int dy = 0; dy < 2; ++dy) {
    for(int dx = 0; dx < 2; ++dx) {
      frame.set_sample(static_cast<std::size_t>(star_x + dx),
                       static_cast<std::size_t>(star_y + dy),
                       astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = false});
    }
  }
  return frame;
}

void estimates_integer_translation_from_cfa_safe_stars() {
  const astrocfa::CfaFrame reference = synthetic_star_field(16, 16, 6, 6);
  const astrocfa::CfaFrame moving = synthetic_star_field(16, 16, 8, 4);

  const astrocfa::RegistrationResult result = astrocfa::estimate_integer_star_translation(
      reference, moving,
      astrocfa::RegistrationOptions{.max_shift_pixels = 8, .sigma_threshold = 2.0, .sample_stride = 1});

  require(result.offset.dx == -2.0, "Moving frame should shift left by two CFA pixels");
  require(result.offset.dy == 2.0, "Moving frame should shift down by two CFA pixels");
  require(result.score > 0.0, "Registration score should be positive");
}

void returns_zero_when_no_stellar_signal_matches() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame reference(16, 16, rggb);
  astrocfa::CfaFrame moving(16, 16, rggb);
  for(std::size_t y = 0; y < reference.height(); ++y) {
    for(std::size_t x = 0; x < reference.width(); ++x) {
      reference.set_sample(
          x, y, astrocfa::CfaSample{.value = 0.01F, .valid = true, .clipped = false});
      moving.set_sample(
          x, y, astrocfa::CfaSample{.value = 0.01F, .valid = true, .clipped = false});
    }
  }

  const astrocfa::RegistrationResult result = astrocfa::estimate_integer_star_translation(
      reference, moving,
      astrocfa::RegistrationOptions{.max_shift_pixels = 8, .sigma_threshold = 2.0, .sample_stride = 1});

  require(result.offset.dx == 0.0, "No stellar signal should keep zero dx");
  require(result.offset.dy == 0.0, "No stellar signal should keep zero dy");
  require(result.score == 0.0, "No stellar signal should report zero score");
}

} // namespace

int main() {
  try {
    estimates_integer_translation_from_cfa_safe_stars();
    returns_zero_when_no_stellar_signal_matches();
  } catch(const std::exception &error) {
    std::cerr << "registration_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
