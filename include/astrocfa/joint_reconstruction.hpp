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
  double psf_sigma = 0.0;
};

struct JointReconstructionOptions {
  std::size_t scale = 1;
  std::size_t iterations = 6;
  double learning_rate = 0.85;
  double huber_sigma = 4.0;
  double luma_smoothness = 0.2;
  double chroma_smoothness = 0.9;
  double edge_sensitivity = 24.0;
  bool stop_on_discrepancy = true;
  std::size_t minimum_iterations = 2;
  double discrepancy_target = 1.0;
  double discrepancy_tolerance = 0.1;
  NoiseModel noise;
};

struct JointReconstructionStats {
  std::size_t frames = 0;
  std::size_t measurements = 0;
  std::size_t iterations = 0;
  std::size_t maximum_iterations = 0;
  std::size_t robust_outliers = 0;
  std::size_t psf_frames = 0;
  double minimum_psf_sigma = 0.0;
  double maximum_psf_sigma = 0.0;
  double initial_rmse = 0.0;
  double final_rmse = 0.0;
  double final_normalized_mae = 0.0;
  double initial_reduced_chi_square = 0.0;
  double final_reduced_chi_square = 0.0;
  bool stopped_by_discrepancy = false;
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
