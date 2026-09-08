#pragma once

#include "astrocfa/cfa.hpp"

#include <string>

namespace astrocfa {

enum class RgbTransfer {
  linear,
  srgb,
};

[[nodiscard]] RgbImage read_rgb_tiff(const std::string &path,
                                     RgbTransfer transfer = RgbTransfer::linear);

} // namespace astrocfa
