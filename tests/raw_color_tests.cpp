#include "astrocfa/raw_color.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

bool near(double left, double right) {
  return std::abs(left - right) < 1.0e-6;
}

astrocfa::RawColorMetadata test_metadata() {
  return astrocfa::RawColorMetadata{
      .has_as_shot_white_balance = true,
      .has_daylight_white_balance = true,
      .has_camera_to_srgb = true,
      .as_shot_white_balance = {2.0, 1.0, 1.5},
      .daylight_white_balance = {1.8, 1.0, 1.3},
      .camera_to_srgb = {
          std::array<double, 3>{1.0, 0.1, 0.0},
          std::array<double, 3>{0.0, 0.9, 0.0},
          std::array<double, 3>{0.0, 0.0, 1.1},
      },
  };
}

void applies_white_balance_before_camera_matrix() {
  astrocfa::RgbImage input(1, 1);
  input.set_pixel(0, 0, {.r = 0.1F, .g = 0.2F, .b = 0.3F});
  const astrocfa::RawColorResult result =
      astrocfa::apply_raw_color(input, test_metadata());
  const astrocfa::RgbPixel pixel = result.image.pixel(0, 0);
  require(near(pixel.r, 0.22), "Red should use balanced camera RGB");
  require(near(pixel.g, 0.18), "Green should use the camera matrix");
  require(near(pixel.b, 0.495), "Blue should use balanced camera RGB");
  require(result.stats.used_metadata_white_balance,
          "As-shot mode should report metadata use");
  require(result.stats.converted_to_srgb,
          "Camera matrix application should be reported");
}

void unity_camera_space_preserves_input() {
  astrocfa::RgbImage input(1, 1);
  input.set_pixel(0, 0, {.r = 0.1F, .g = 0.2F, .b = 0.3F});
  const astrocfa::RawColorResult result = astrocfa::apply_raw_color(
      input, {}, astrocfa::RawColorOptions{
                     .white_balance = astrocfa::WhiteBalanceMode::unity,
                     .convert_to_srgb = false,
                 });
  const astrocfa::RgbPixel pixel = result.image.pixel(0, 0);
  require(near(pixel.r, 0.1) && near(pixel.g, 0.2) && near(pixel.b, 0.3),
          "Unity camera-space development should preserve sensor RGB");
}

void custom_balance_is_green_normalized_and_unclipped() {
  astrocfa::RgbImage input(1, 1);
  input.set_pixel(0, 0, {.r = 0.8F, .g = 0.2F, .b = 0.4F});
  const astrocfa::RawColorResult result = astrocfa::apply_raw_color(
      input, {}, astrocfa::RawColorOptions{
                     .white_balance = astrocfa::WhiteBalanceMode::custom,
                     .custom_white_balance = {4.0, 2.0, 1.0},
                     .convert_to_srgb = false,
                 });
  const astrocfa::RgbPixel pixel = result.image.pixel(0, 0);
  require(near(pixel.r, 1.6) && near(pixel.g, 0.2) && near(pixel.b, 0.2),
          "Custom WB should normalize multipliers around green");
  require(result.stats.over_range_pixels == 1,
          "Linear development should report, not clip, values over one");
}

void missing_requested_metadata_is_rejected() {
  bool rejected = false;
  try {
    (void)astrocfa::apply_raw_color(astrocfa::RgbImage(1, 1), {});
  } catch(const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "Missing as-shot WB should fail explicitly");
}

void automatic_mode_uses_available_metadata_conservatively() {
  astrocfa::RawColorMetadata metadata;
  metadata.has_daylight_white_balance = true;
  const astrocfa::RawColorOptions options =
      astrocfa::automatic_raw_color_options(metadata);
  require(options.white_balance == astrocfa::WhiteBalanceMode::daylight,
          "Auto color should fall back from missing as-shot WB to daylight");
  require(!options.convert_to_srgb,
          "Auto color should preserve camera RGB without a valid matrix");
}

} // namespace

int main() {
  try {
    applies_white_balance_before_camera_matrix();
    unity_camera_space_preserves_input();
    custom_balance_is_green_normalized_and_unclipped();
    missing_requested_metadata_is_rejected();
    automatic_mode_uses_available_metadata_conservatively();
  } catch(const std::exception &error) {
    std::cerr << "raw_color_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
