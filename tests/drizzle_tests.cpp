#include "astrocfa/drizzle.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void require_near(double actual, double expected, double tolerance, const char *message) {
  if(std::abs(actual - expected) > tolerance) {
    std::cerr << message << ": actual=" << actual << " expected=" << expected << "\n";
    throw std::runtime_error(message);
  }
}

void single_sample_deposits_to_matching_phase() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(2, 2, rggb);
  frame.set_sample(0, 0, astrocfa::CfaSample{.value = 0.8F, .valid = true, .clipped = false});
  frame.set_sample(1, 0, astrocfa::CfaSample{.value = 0.2F, .valid = false, .clipped = false});
  frame.set_sample(0, 1, astrocfa::CfaSample{.value = 0.2F, .valid = false, .clipped = false});
  frame.set_sample(1, 1, astrocfa::CfaSample{.value = 0.2F, .valid = false, .clipped = false});

  astrocfa::CfaDrizzleAccumulator drizzle(2, 2, rggb,
                                          astrocfa::DrizzleOptions{.scale = 2, .drop_shrink = 1.0});
  drizzle.add_frame(frame, astrocfa::SubpixelOffset{});

  const astrocfa::DrizzleCoverageStats red = drizzle.coverage(astrocfa::CfaColor::red);
  const astrocfa::DrizzleCoverageStats blue = drizzle.coverage(astrocfa::CfaColor::blue);
  require(red.covered_samples > 0, "Red phase should receive the valid red sample");
  require(blue.covered_samples == 0, "Blue phase should not receive red samples");
}

void matching_samples_average_by_weight() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame first(2, 2, rggb);
  astrocfa::CfaFrame second(2, 2, rggb);
  first.set_sample(0, 0, astrocfa::CfaSample{.value = 0.25F, .valid = true, .clipped = false});
  second.set_sample(0, 0, astrocfa::CfaSample{.value = 0.75F, .valid = true, .clipped = false});
  for(std::size_t y = 0; y < 2; ++y) {
    for(std::size_t x = 0; x < 2; ++x) {
      if(x == 0 && y == 0) {
        continue;
      }
      first.set_valid(x, y, false);
      second.set_valid(x, y, false);
    }
  }

  astrocfa::CfaDrizzleAccumulator drizzle(2, 2, rggb,
                                          astrocfa::DrizzleOptions{.scale = 2, .drop_shrink = 1.0});
  drizzle.add_frame(first, astrocfa::SubpixelOffset{}, 1.0);
  drizzle.add_frame(second, astrocfa::SubpixelOffset{}, 3.0);

  bool found = false;
  for(std::size_t y = 0; y < drizzle.output_height(); ++y) {
    for(std::size_t x = 0; x < drizzle.output_width(); ++x) {
      const astrocfa::DrizzlePlaneSample sample = drizzle.plane_sample(astrocfa::CfaColor::red, x, y);
      if(sample.weight > 0.0) {
        require_near(sample.value, 0.625, 1.0e-6, "Weighted drizzle average");
        found = true;
      }
    }
  }
  require(found, "At least one drizzle output sample should be covered");
}

void clipped_samples_are_not_deposited() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(2, 2, rggb);
  frame.set_sample(0, 0, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});

  astrocfa::CfaDrizzleAccumulator drizzle(2, 2, rggb,
                                          astrocfa::DrizzleOptions{.scale = 2, .drop_shrink = 1.0});
  drizzle.add_frame(frame, astrocfa::SubpixelOffset{});
  require(drizzle.coverage(astrocfa::CfaColor::red).covered_samples == 0,
          "Clipped samples should not be accumulated");
}

} // namespace

int main() {
  try {
    single_sample_deposits_to_matching_phase();
    matching_samples_average_by_weight();
    clipped_samples_are_not_deposited();
  } catch(const std::exception &error) {
    std::cerr << "drizzle_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}

