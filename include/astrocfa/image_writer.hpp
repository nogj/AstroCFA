#pragma once

#include "astrocfa/cfa.hpp"

#include <string>

namespace astrocfa {

struct ImageWriteOptions {
  int jpeg_quality = 92;
  bool apply_preview_gamma = true;
};

void write_rgb_image(const RgbImage &image, const std::string &path,
                     ImageWriteOptions options = {});
void write_rgb_tiff16(const RgbImage &image, const std::string &path);
void write_linear_cfa_dng16(const CfaFrame &cfa, const std::string &path);
void write_rgb_jpeg8(const RgbImage &image, const std::string &path,
                     ImageWriteOptions options = {});

} // namespace astrocfa
