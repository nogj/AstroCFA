#pragma once

#include "astrocfa/cfa.hpp"

#include <array>
#include <cstddef>

namespace astrocfa {

using ColorMatrix3x3 = std::array<std::array<double, 3>, 3>;

struct RawColorMetadata {
  bool has_as_shot_white_balance = false;
  bool has_daylight_white_balance = false;
  bool has_camera_to_srgb = false;
  std::array<double, 3> as_shot_white_balance = {1.0, 1.0, 1.0};
  std::array<double, 3> daylight_white_balance = {1.0, 1.0, 1.0};
  ColorMatrix3x3 camera_to_srgb = {
      std::array<double, 3>{1.0, 0.0, 0.0},
      std::array<double, 3>{0.0, 1.0, 0.0},
      std::array<double, 3>{0.0, 0.0, 1.0},
  };
};

enum class WhiteBalanceMode {
  as_shot,
  daylight,
  unity,
  custom,
};

struct RawColorOptions {
  WhiteBalanceMode white_balance = WhiteBalanceMode::as_shot;
  std::array<double, 3> custom_white_balance = {1.0, 1.0, 1.0};
  bool convert_to_srgb = true;
};

struct RawColorStats {
  std::array<double, 3> white_balance = {1.0, 1.0, 1.0};
  bool used_metadata_white_balance = false;
  bool converted_to_srgb = false;
  std::size_t negative_pixels = 0;
  std::size_t over_range_pixels = 0;
  double maximum_component = 0.0;
};

struct RawColorResult {
  RgbImage image;
  RawColorStats stats;
};

[[nodiscard]] RawColorResult apply_raw_color(
    const RgbImage &camera_rgb, const RawColorMetadata &metadata,
    RawColorOptions options = {});
[[nodiscard]] RawColorOptions automatic_raw_color_options(
    const RawColorMetadata &metadata);

} // namespace astrocfa
