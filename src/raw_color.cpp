#include "astrocfa/raw_color.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

std::array<double, 3> normalized_white_balance(
    std::array<double, 3> multipliers) {
  for(double multiplier : multipliers) {
    if(!std::isfinite(multiplier) || multiplier <= 0.0) {
      throw std::invalid_argument(
          "White-balance multipliers must be finite and positive");
    }
  }
  const double green = multipliers[1];
  for(double &multiplier : multipliers) {
    multiplier /= green;
  }
  return multipliers;
}

bool finite_matrix(const astrocfa::ColorMatrix3x3 &matrix) {
  for(const auto &row : matrix) {
    for(double value : row) {
      if(!std::isfinite(value)) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

namespace astrocfa {

RawColorOptions automatic_raw_color_options(const RawColorMetadata &metadata) {
  RawColorOptions options;
  if(metadata.has_as_shot_white_balance) {
    options.white_balance = WhiteBalanceMode::as_shot;
  } else if(metadata.has_daylight_white_balance) {
    options.white_balance = WhiteBalanceMode::daylight;
  } else {
    options.white_balance = WhiteBalanceMode::unity;
  }
  options.convert_to_srgb = metadata.has_camera_to_srgb;
  return options;
}

RawColorResult apply_raw_color(const RgbImage &camera_rgb,
                               const RawColorMetadata &metadata,
                               RawColorOptions options) {
  std::array<double, 3> white_balance = {1.0, 1.0, 1.0};
  bool used_metadata_white_balance = false;
  switch(options.white_balance) {
  case WhiteBalanceMode::as_shot:
    if(!metadata.has_as_shot_white_balance) {
      throw std::invalid_argument("RAW has no as-shot white balance metadata");
    }
    white_balance = metadata.as_shot_white_balance;
    used_metadata_white_balance = true;
    break;
  case WhiteBalanceMode::daylight:
    if(!metadata.has_daylight_white_balance) {
      throw std::invalid_argument("RAW has no daylight white balance metadata");
    }
    white_balance = metadata.daylight_white_balance;
    used_metadata_white_balance = true;
    break;
  case WhiteBalanceMode::unity:
    break;
  case WhiteBalanceMode::custom:
    white_balance = options.custom_white_balance;
    break;
  }
  white_balance = normalized_white_balance(white_balance);

  if(options.convert_to_srgb &&
     (!metadata.has_camera_to_srgb || !finite_matrix(metadata.camera_to_srgb))) {
    throw std::invalid_argument("RAW has no usable camera-to-sRGB matrix");
  }

  RawColorResult result{
      .image = RgbImage(camera_rgb.width(), camera_rgb.height()),
      .stats = RawColorStats{
          .white_balance = white_balance,
          .used_metadata_white_balance = used_metadata_white_balance,
          .converted_to_srgb = options.convert_to_srgb,
      },
  };
  for(std::size_t y = 0; y < camera_rgb.height(); ++y) {
    for(std::size_t x = 0; x < camera_rgb.width(); ++x) {
      const RgbPixel source = camera_rgb.pixel(x, y);
      const std::array<double, 3> balanced = {
          source.r * white_balance[0],
          source.g * white_balance[1],
          source.b * white_balance[2],
      };
      std::array<double, 3> developed = balanced;
      if(options.convert_to_srgb) {
        for(std::size_t output = 0; output < 3; ++output) {
          developed[output] = 0.0;
          for(std::size_t input = 0; input < 3; ++input) {
            developed[output] +=
                metadata.camera_to_srgb[output][input] * balanced[input];
          }
        }
      }
      result.image.set_pixel(
          x, y,
          RgbPixel{.r = static_cast<float>(developed[0]),
                   .g = static_cast<float>(developed[1]),
                   .b = static_cast<float>(developed[2])});
      const double minimum =
          std::min({developed[0], developed[1], developed[2]});
      const double maximum =
          std::max({developed[0], developed[1], developed[2]});
      result.stats.negative_pixels += minimum < 0.0 ? 1U : 0U;
      result.stats.over_range_pixels += maximum > 1.0 ? 1U : 0U;
      result.stats.maximum_component =
          std::max(result.stats.maximum_component, maximum);
    }
  }
  return result;
}

} // namespace astrocfa
