#include "astrocfa/image_reader.hpp"

#include <tiffio.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template <typename Sample>
float normalized_integer(Sample value) {
  return static_cast<float>(value) /
         static_cast<float>(std::numeric_limits<Sample>::max());
}

} // namespace

namespace astrocfa {

RgbImage read_linear_rgb_tiff(const std::string &path) {
  TIFF *tiff = TIFFOpen(path.c_str(), "r");
  if(tiff == nullptr) {
    throw std::runtime_error("Cannot open candidate TIFF: " + path);
  }

  std::uint32_t width = 0;
  std::uint32_t height = 0;
  std::uint16_t channels = 0;
  std::uint16_t bits = 0;
  std::uint16_t planar = PLANARCONFIG_CONTIG;
  std::uint16_t photometric = 0;
  std::uint16_t sample_format = SAMPLEFORMAT_UINT;
  std::uint16_t orientation = ORIENTATION_TOPLEFT;
  TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
  TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &channels);
  TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits);
  TIFFGetFieldDefaulted(tiff, TIFFTAG_PLANARCONFIG, &planar);
  TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
  TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLEFORMAT, &sample_format);
  TIFFGetFieldDefaulted(tiff, TIFFTAG_ORIENTATION, &orientation);

  const bool supported_integer = sample_format == SAMPLEFORMAT_UINT &&
                                 (bits == 8 || bits == 16);
  const bool supported_float = sample_format == SAMPLEFORMAT_IEEEFP && bits == 32;
  if(width == 0 || height == 0 || channels < 3 || photometric != PHOTOMETRIC_RGB ||
     planar != PLANARCONFIG_CONTIG || orientation != ORIENTATION_TOPLEFT ||
     (!supported_integer && !supported_float)) {
    TIFFClose(tiff);
    throw std::invalid_argument(
        "Candidate must be a top-left, contiguous linear RGB TIFF with unsigned "
        "8/16-bit or float32 samples");
  }

  RgbImage image(width, height);
  std::vector<std::uint8_t> row(static_cast<std::size_t>(TIFFScanlineSize(tiff)));
  for(std::uint32_t y = 0; y < height; ++y) {
    if(TIFFReadScanline(tiff, row.data(), y, 0) < 0) {
      TIFFClose(tiff);
      throw std::runtime_error("Failed while reading candidate TIFF scanline");
    }
    for(std::uint32_t x = 0; x < width; ++x) {
      const std::size_t offset = static_cast<std::size_t>(x) * channels;
      RgbPixel pixel;
      if(bits == 8) {
        const auto *samples = reinterpret_cast<const std::uint8_t *>(row.data());
        pixel = {.r = normalized_integer(samples[offset]),
                 .g = normalized_integer(samples[offset + 1U]),
                 .b = normalized_integer(samples[offset + 2U])};
      } else if(bits == 16) {
        const auto *samples = reinterpret_cast<const std::uint16_t *>(row.data());
        pixel = {.r = normalized_integer(samples[offset]),
                 .g = normalized_integer(samples[offset + 1U]),
                 .b = normalized_integer(samples[offset + 2U])};
      } else {
        const auto *samples = reinterpret_cast<const float *>(row.data());
        pixel = {.r = samples[offset], .g = samples[offset + 1U],
                 .b = samples[offset + 2U]};
      }
      if(!std::isfinite(pixel.r) || !std::isfinite(pixel.g) ||
         !std::isfinite(pixel.b) || pixel.r < 0.0F || pixel.r > 1.0F ||
         pixel.g < 0.0F || pixel.g > 1.0F || pixel.b < 0.0F || pixel.b > 1.0F) {
        TIFFClose(tiff);
        throw std::invalid_argument(
            "Candidate TIFF contains non-finite or out-of-range linear RGB samples");
      }
      image.set_pixel(x, y, pixel);
    }
  }

  TIFFClose(tiff);
  return image;
}

} // namespace astrocfa
