#include "astrocfa/image_writer.hpp"

#include <jpeglib.h>
#include <tiffio.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

float clamp_display(float value) {
  return std::clamp(value, 0.0F, 1.0F);
}

std::uint16_t to_u16(float value) {
  return static_cast<std::uint16_t>(std::lround(clamp_display(value) * 65535.0F));
}

unsigned char to_u8(float value, bool gamma) {
  float mapped = clamp_display(value);
  if(gamma) {
    mapped = std::pow(mapped, 1.0F / 2.2F);
  }
  return static_cast<unsigned char>(std::lround(mapped * 255.0F));
}

std::string lowercase_extension(const std::string &path) {
  const std::size_t dot = path.find_last_of('.');
  if(dot == std::string::npos) {
    return {};
  }

  std::string extension = path.substr(dot + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
  return extension;
}

} // namespace

namespace astrocfa {

void write_rgb_image(const RgbImage &image, const std::string &path,
                     ImageWriteOptions options) {
  const std::string extension = lowercase_extension(path);
  if(extension == "tif" || extension == "tiff") {
    write_rgb_tiff16(image, path);
    return;
  }
  if(extension == "jpg" || extension == "jpeg") {
    write_rgb_jpeg8(image, path, options);
    return;
  }

  throw std::invalid_argument("Unsupported output extension: " + extension);
}

void write_rgb_tiff16(const RgbImage &image, const std::string &path) {
  TIFF *tiff = TIFFOpen(path.c_str(), "w");
  if(tiff == nullptr) {
    throw std::runtime_error("Cannot open TIFF for writing: " + path);
  }

  TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, static_cast<std::uint32_t>(image.width()));
  TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, static_cast<std::uint32_t>(image.height()));
  TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
  TIFFSetField(tiff, TIFFTAG_PREDICTOR, PREDICTOR_HORIZONTAL);
  TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tiff, 0));

  std::vector<std::uint16_t> row(image.width() * 3U);
  for(std::size_t y = 0; y < image.height(); ++y) {
    for(std::size_t x = 0; x < image.width(); ++x) {
      const RgbPixel pixel = image.pixel(x, y);
      row[x * 3U + 0U] = to_u16(pixel.r);
      row[x * 3U + 1U] = to_u16(pixel.g);
      row[x * 3U + 2U] = to_u16(pixel.b);
    }

    if(TIFFWriteScanline(tiff, row.data(), static_cast<std::uint32_t>(y), 0) < 0) {
      TIFFClose(tiff);
      throw std::runtime_error("Failed while writing TIFF scanline");
    }
  }

  TIFFClose(tiff);
}

void write_linear_cfa_dng16(const CfaFrame &cfa, const std::string &path) {
  TIFF *tiff = TIFFOpen(path.c_str(), "w");
  if(tiff == nullptr) {
    throw std::runtime_error("Cannot open DNG for writing: " + path);
  }

  const std::array<std::uint8_t, 4> dng_version = {1, 4, 0, 0};
  const std::array<std::uint8_t, 4> backward_version = {1, 1, 0, 0};
  const std::array<std::uint16_t, 2> cfa_dimensions = {2, 2};
  std::array<std::uint8_t, 4> cfa_pattern{};
  for(std::size_t index = 0; index < cfa_pattern.size(); ++index) {
    const CfaColor color = cfa.pattern().phase[index];
    cfa_pattern[index] = color == CfaColor::red ? 0U
                         : color == CfaColor::blue ? 2U
                                                   : 1U;
  }
  const std::array<std::uint8_t, 3> cfa_plane_color = {0, 1, 2};
  const std::array<float, 9> color_matrix = {
      1.0F, 0.0F, 0.0F,
      0.0F, 1.0F, 0.0F,
      0.0F, 0.0F, 1.0F,
  };
  const std::array<float, 3> as_shot_neutral = {1.0F, 1.0F, 1.0F};
  const std::array<std::uint32_t, 1> white_level = {65535U};

  TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, static_cast<std::uint32_t>(cfa.width()));
  TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, static_cast<std::uint32_t>(cfa.height()));
  TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
  TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
  TIFFSetField(tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_CFA);
  TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tiff, 0));
  TIFFSetField(tiff, TIFFTAG_MAKE, "AstroCFA");
  TIFFSetField(tiff, TIFFTAG_MODEL, "Synthetic CFA benchmark");
  TIFFSetField(tiff, TIFFTAG_SOFTWARE, "AstroCFA benchmark-debayer");
  TIFFSetField(tiff, TIFFTAG_DNGVERSION, dng_version.data());
  TIFFSetField(tiff, TIFFTAG_DNGBACKWARDVERSION, backward_version.data());
  TIFFSetField(tiff, TIFFTAG_UNIQUECAMERAMODEL, "AstroCFA Synthetic Linear CFA");
  TIFFSetField(tiff, TIFFTAG_CFAREPEATPATTERNDIM, cfa_dimensions.data());
  TIFFSetField(tiff, TIFFTAG_CFAPATTERN,
               static_cast<std::uint16_t>(cfa_pattern.size()), cfa_pattern.data());
  TIFFSetField(tiff, TIFFTAG_CFAPLANECOLOR,
               static_cast<std::uint16_t>(cfa_plane_color.size()),
               cfa_plane_color.data());
  TIFFSetField(tiff, TIFFTAG_CFALAYOUT, 1);
  TIFFSetField(tiff, TIFFTAG_WHITELEVEL, 1U, white_level.data());
  TIFFSetField(tiff, TIFFTAG_COLORMATRIX1,
               static_cast<std::uint16_t>(color_matrix.size()), color_matrix.data());
  TIFFSetField(tiff, TIFFTAG_ASSHOTNEUTRAL,
               static_cast<std::uint16_t>(as_shot_neutral.size()),
               as_shot_neutral.data());
  TIFFSetField(tiff, TIFFTAG_CALIBRATIONILLUMINANT1, 21);

  std::vector<std::uint16_t> row(cfa.width());
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      row[x] = to_u16(cfa.sample_info(x, y).value);
    }
    if(TIFFWriteScanline(tiff, row.data(), static_cast<std::uint32_t>(y), 0) < 0) {
      TIFFClose(tiff);
      throw std::runtime_error("Failed while writing DNG scanline");
    }
  }

  TIFFClose(tiff);
}

void write_rgb_jpeg8(const RgbImage &image, const std::string &path,
                     ImageWriteOptions options) {
  FILE *file = std::fopen(path.c_str(), "wb");
  if(file == nullptr) {
    throw std::runtime_error("Cannot open JPEG for writing: " + path);
  }

  jpeg_compress_struct compressor{};
  jpeg_error_mgr error_manager{};
  compressor.err = jpeg_std_error(&error_manager);
  jpeg_create_compress(&compressor);
  jpeg_stdio_dest(&compressor, file);

  compressor.image_width = static_cast<JDIMENSION>(image.width());
  compressor.image_height = static_cast<JDIMENSION>(image.height());
  compressor.input_components = 3;
  compressor.in_color_space = JCS_RGB;
  jpeg_set_defaults(&compressor);
  jpeg_set_quality(&compressor, std::clamp(options.jpeg_quality, 1, 100), TRUE);
  jpeg_start_compress(&compressor, TRUE);

  std::vector<unsigned char> row(image.width() * 3U);
  while(compressor.next_scanline < compressor.image_height) {
    const std::size_t y = compressor.next_scanline;
    for(std::size_t x = 0; x < image.width(); ++x) {
      const RgbPixel pixel = image.pixel(x, y);
      row[x * 3U + 0U] = to_u8(pixel.r, options.apply_preview_gamma);
      row[x * 3U + 1U] = to_u8(pixel.g, options.apply_preview_gamma);
      row[x * 3U + 2U] = to_u8(pixel.b, options.apply_preview_gamma);
    }

    JSAMPROW row_pointer = row.data();
    jpeg_write_scanlines(&compressor, &row_pointer, 1);
  }

  jpeg_finish_compress(&compressor);
  jpeg_destroy_compress(&compressor);
  std::fclose(file);
}

} // namespace astrocfa
