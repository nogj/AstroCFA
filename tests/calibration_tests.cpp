#include "astrocfa/calibration.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

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

astrocfa::CfaFrame filled_frame(float value) {
  astrocfa::CfaFrame frame(4, 4, astrocfa::BayerPattern{});
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{
                                 .value = value,
                                 .valid = true,
                                 .clipped = false,
                             });
    }
  }
  return frame;
}

astrocfa::CfaFrame filled_frame(std::size_t width, std::size_t height, float value) {
  astrocfa::CfaFrame frame(width, height, astrocfa::BayerPattern{});
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true});
    }
  }
  return frame;
}

void subtracts_dark_without_double_subtracting_bias_by_default() {
  astrocfa::CfaFrame light = filled_frame(0.50F);
  astrocfa::CfaFrame bias = filled_frame(0.03F);
  astrocfa::CfaFrame dark = filled_frame(0.10F);

  const astrocfa::CalibrationResult result =
      astrocfa::calibrate_cfa(light, astrocfa::CalibrationInputs{
                                         .bias = &bias,
                                         .dark = &dark,
                                     });
  require_near(result.cfa.sample(0, 0), 0.40, 1.0e-6,
               "Master dark should include bias by default");
}

void can_subtract_bias_when_dark_is_bias_free() {
  astrocfa::CfaFrame light = filled_frame(0.50F);
  astrocfa::CfaFrame bias = filled_frame(0.03F);
  astrocfa::CfaFrame dark = filled_frame(0.10F);

  const astrocfa::CalibrationResult result =
      astrocfa::calibrate_cfa(light, astrocfa::CalibrationInputs{
                                         .bias = &bias,
                                         .dark = &dark,
                                         .options = {.dark_includes_bias = false},
                                     });
  require_near(result.cfa.sample(0, 0), 0.37, 1.0e-6,
               "Bias-free dark should allow explicit bias subtraction");
}

void normalizes_flat_per_cfa_phase() {
  astrocfa::CfaFrame light = filled_frame(0.50F);
  astrocfa::CfaFrame flat(4, 4, astrocfa::BayerPattern{});
  for(std::size_t y = 0; y < flat.height(); ++y) {
    for(std::size_t x = 0; x < flat.width(); ++x) {
      const float value = x < 2 ? 0.25F : 0.50F;
      flat.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true});
    }
  }

  const astrocfa::CalibrationResult result =
      astrocfa::calibrate_cfa(light, astrocfa::CalibrationInputs{.flat = &flat});
  require(result.stats.flat_floor_samples == 0, "Valid flat should avoid floor");
  require(result.cfa.sample(0, 0) > result.cfa.sample(2, 0),
          "Darker flat region should brighten calibrated light");
}

void builds_median_master_and_rejects_outlier() {
  std::vector<astrocfa::CfaFrame> frames;
  frames.push_back(filled_frame(0.10F));
  frames.push_back(filled_frame(0.11F));
  frames.push_back(filled_frame(0.12F));
  frames.push_back(filled_frame(0.13F));
  frames.push_back(filled_frame(0.95F));

  const astrocfa::MasterBuildResult master =
      astrocfa::build_master_cfa(frames, astrocfa::MasterBuildOptions{
                                             .sigma_clip = 2.0,
                                             .min_clip_samples = 5,
                                         });
  require(master.stats.frames == 5, "Master stats should count frames");
  require(master.stats.rejected_samples > 0, "Outlier samples should be rejected");
  require_near(master.cfa.sample(0, 0), 0.115, 0.01,
               "Robust master should ignore strong outlier");
}

void master_builder_rejects_incompatible_dimensions() {
  std::vector<astrocfa::CfaFrame> frames;
  frames.push_back(filled_frame(0.10F));
  frames.emplace_back(2, 2, astrocfa::BayerPattern{});

  bool threw = false;
  try {
    (void)astrocfa::build_master_cfa(frames);
  } catch(const std::invalid_argument &) {
    threw = true;
  }
  require(threw, "Master builder should reject incompatible frame dimensions");
}

void detects_and_repairs_master_supported_sensor_defects() {
  astrocfa::CfaFrame light = filled_frame(12, 12, 0.50F);
  astrocfa::CfaFrame dark = filled_frame(12, 12, 0.01F);
  astrocfa::CfaFrame flat = filled_frame(12, 12, 0.50F);
  dark.set_sample(6, 6, 0.40F);
  flat.set_sample(7, 7, 0.02F);

  const astrocfa::CalibrationResult result = astrocfa::calibrate_cfa(
      light, astrocfa::CalibrationInputs{.dark = &dark, .flat = &flat});

  require(result.defects.has(6, 6, astrocfa::SensorDefect::hot),
          "Dark outlier should be marked hot");
  require(result.defects.has(7, 7, astrocfa::SensorDefect::dead),
          "Low flat response should be marked dead");
  require(result.stats.hot_pixels == 1, "Hot pixel stats should be exact");
  require(result.stats.dead_pixels == 1, "Dead pixel stats should be exact");
  require(result.stats.repaired_pixels == 2, "Both defects should be repaired");
  require(result.stats.unrepaired_pixels == 0, "Interior defects should have neighbors");
  require_near(result.cfa.sample(6, 6), 0.49, 1.0e-5,
               "Hot pixel should use same-phase calibrated neighbors");
  require_near(result.cfa.sample(7, 7), 0.49, 1.0e-5,
               "Dead pixel should use same-phase calibrated neighbors");
}

void can_disable_cosmetic_correction() {
  astrocfa::CfaFrame light = filled_frame(12, 12, 0.50F);
  astrocfa::CfaFrame dark = filled_frame(12, 12, 0.01F);
  dark.set_sample(6, 6, 0.40F);

  const astrocfa::CalibrationResult result = astrocfa::calibrate_cfa(
      light, astrocfa::CalibrationInputs{
                 .dark = &dark,
                 .options = {.cosmetic = {.enabled = false}},
             });
  require(!result.defects.defective(6, 6), "Disabled correction should skip detection");
  require_near(result.cfa.sample(6, 6), 0.10, 1.0e-5,
               "Disabled correction should preserve ordinary calibration result");
}

} // namespace

int main() {
  try {
    subtracts_dark_without_double_subtracting_bias_by_default();
    can_subtract_bias_when_dark_is_bias_free();
    normalizes_flat_per_cfa_phase();
    builds_median_master_and_rejects_outlier();
    master_builder_rejects_incompatible_dimensions();
    detects_and_repairs_master_supported_sensor_defects();
    can_disable_cosmetic_correction();
  } catch(const std::exception &error) {
    std::cerr << "calibration_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
