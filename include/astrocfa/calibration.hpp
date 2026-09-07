#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>
#include <vector>

namespace astrocfa {

struct CalibrationOptions {
  bool dark_includes_bias = true;
  double flat_floor = 0.05;
};

struct CalibrationInputs {
  const CfaFrame *bias = nullptr;
  const CfaFrame *dark = nullptr;
  const CfaFrame *flat = nullptr;
  CalibrationOptions options;
};

struct CalibrationStats {
  std::size_t samples = 0;
  std::size_t clipped_samples = 0;
  std::size_t invalid_samples = 0;
  std::size_t flat_floor_samples = 0;
  double mean_before = 0.0;
  double mean_after = 0.0;
  double mean_bias_subtracted = 0.0;
  double mean_dark_subtracted = 0.0;
  double flat_phase_mean[4] = {1.0, 1.0, 1.0, 1.0};
};

struct CalibrationResult {
  CfaFrame cfa;
  CalibrationStats stats;
};

struct MasterBuildOptions {
  double sigma_clip = 4.5;
  std::size_t min_clip_samples = 5;
};

struct MasterBuildStats {
  std::size_t frames = 0;
  std::size_t samples = 0;
  std::size_t invalid_output_samples = 0;
  std::size_t rejected_samples = 0;
  double mean = 0.0;
};

struct MasterBuildResult {
  CfaFrame cfa;
  MasterBuildStats stats;
};

[[nodiscard]] CalibrationResult calibrate_cfa(const CfaFrame &light,
                                              CalibrationInputs inputs = {});
[[nodiscard]] MasterBuildResult build_master_cfa(
    const std::vector<CfaFrame> &frames, MasterBuildOptions options = {});

} // namespace astrocfa
