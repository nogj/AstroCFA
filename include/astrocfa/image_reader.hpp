#pragma once

#include "astrocfa/cfa.hpp"

#include <string>

namespace astrocfa {

[[nodiscard]] RgbImage read_linear_rgb_tiff(const std::string &path);

} // namespace astrocfa
