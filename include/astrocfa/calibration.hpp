#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace astrocfa {

enum class SensorDefect : std::uint8_t {
  none = 0,
  hot = 1U << 0U,
  dead = 1U << 1U,
  invalid_master = 1U << 2U,
};

class DefectMap {
public:
  DefectMap(std::size_t width, std::size_t height);

  [[nodiscard]] std::size_t width() const { return width_; }
  [[nodiscard]] std::size_t height() const { return height_; }
  [[nodiscard]] SensorDefect at(std::size_t x, std::size_t y) const;
  [[nodiscard]] bool has(std::size_t x, std::size_t y, SensorDefect defect) const;
  [[nodiscard]] bool defective(std::size_t x, std::size_t y) const;
  void add(std::size_t x, std::size_t y, SensorDefect defect);

private:
  [[nodiscard]] std::size_t offset(std::size_t x, std::size_t y) const;

  std::size_t width_ = 0;
  std::size_t height_ = 0;
  std::vector<std::uint8_t> flags_;
};

struct CosmeticCorrectionOptions {
  bool enabled = true;
  std::size_t detection_radius = 3;
  std::size_t repair_radius = 3;
  std::size_t min_neighbors = 4;
  double hot_sigma = 8.0;
  double hot_min_excess = 0.01;
  double dead_response_ratio = 0.35;
  double dead_min_deficit = 0.05;
};

struct CalibrationOptions {
  bool dark_includes_bias = true;
  double flat_floor = 0.05;
  CosmeticCorrectionOptions cosmetic;
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
  std::size_t hot_pixels = 0;
  std::size_t dead_pixels = 0;
  std::size_t invalid_master_pixels = 0;
  std::size_t repaired_pixels = 0;
  std::size_t unrepaired_pixels = 0;
  double mean_before = 0.0;
  double mean_after = 0.0;
  double mean_bias_subtracted = 0.0;
  double mean_dark_subtracted = 0.0;
  double flat_phase_mean[4] = {1.0, 1.0, 1.0, 1.0};
};

struct CalibrationResult {
  CfaFrame cfa;
  DefectMap defects;
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
