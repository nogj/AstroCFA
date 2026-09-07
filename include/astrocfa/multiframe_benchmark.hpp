#pragma once

#include "astrocfa/joint_reconstruction.hpp"
#include "astrocfa/reconstruction_metrics.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace astrocfa {

struct MultiframeBenchmarkOptions {
  std::size_t width = 128;
  std::size_t height = 96;
  std::size_t frames = 4;
  std::size_t iterations = 5;
  std::size_t transients_per_frame = 4;
  std::uint32_t seed = 42;
  bool add_noise = true;
  bool vary_seeing = false;
  double luma_smoothness = 0.2;
  double chroma_smoothness = 0.9;
};

struct MultiframeBenchmarkMethod {
  std::string name;
  RgbImage image;
  std::unique_ptr<RgbImage> confidence;
  ReconstructionMetrics metrics;
  JointReconstructionStats solver_stats;
  bool has_solver_stats = false;
};

struct MultiframeBenchmarkResult {
  RgbImage truth;
  std::vector<SyntheticStar> stars;
  std::vector<SubpixelOffset> offsets;
  std::size_t injected_transients = 0;
  std::vector<MultiframeBenchmarkMethod> methods;
};

[[nodiscard]] MultiframeBenchmarkResult run_multiframe_benchmark(
    MultiframeBenchmarkOptions options = {});

} // namespace astrocfa
