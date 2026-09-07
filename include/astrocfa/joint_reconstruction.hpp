#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/drizzle.hpp"
#include "astrocfa/noise_model.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace astrocfa {

struct JointCfaFrame {
  const CfaFrame *cfa = nullptr;
  SubpixelOffset offset;
  double weight = 1.0;
};

struct JointReconstructionOptions {
  std::size_t scale = 1;
  std::size_t iterations = 6;
  double learning_rate = 0.85;
  double huber_sigma = 4.0;
  double chroma_smoothness = 0.08;
  double edge_sensitivity = 24.0;
  NoiseModel noise;
};

struct JointReconstructionStats {
  std::size_t frames = 0;
  std::size_t measurements = 0;
  std::size_t iterations = 0;
  std::size_t robust_outliers = 0;
  double initial_rmse = 0.0;
  double final_rmse = 0.0;
  double final_normalized_mae = 0.0;
  std::array<double, 3> channel_coverage = {0.0, 0.0, 0.0};
};

struct JointReconstructionResult {
  RgbImage image;
  RgbImage confidence;
  JointReconstructionStats stats;
};

[[nodiscard]] JointReconstructionResult reconstruct_joint_cfa(
    const std::vector<JointCfaFrame> &frames,
    JointReconstructionOptions options = {});

} // namespace astrocfa
