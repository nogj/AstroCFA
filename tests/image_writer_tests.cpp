#include "astrocfa/image_reader.hpp"
#include "astrocfa/image_writer.hpp"
#include "astrocfa/raw_loader.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

astrocfa::RgbImage tiny_image() {
  astrocfa::RgbImage image(2, 2);
  image.set_pixel(0, 0, {.r = 0.0F, .g = 0.1F, .b = 0.2F});
  image.set_pixel(1, 0, {.r = 0.3F, .g = 0.4F, .b = 0.5F});
  image.set_pixel(0, 1, {.r = 0.6F, .g = 0.7F, .b = 0.8F});
  image.set_pixel(1, 1, {.r = 0.9F, .g = 1.0F, .b = 1.1F});
  return image;
}

void writes_tiff_and_jpeg() {
  const astrocfa::RgbImage image = tiny_image();
  const std::filesystem::path directory = std::filesystem::temp_directory_path();
  const std::filesystem::path tiff = directory / "astrocfa-image-writer-test.tif";
  const std::filesystem::path jpeg = directory / "astrocfa-image-writer-test.jpg";

  astrocfa::write_rgb_image(image, tiff.string());
  astrocfa::write_rgb_image(image, jpeg.string());

  require(std::filesystem::exists(tiff), "TIFF output should exist");
  require(std::filesystem::file_size(tiff) > 0, "TIFF output should not be empty");
  require(std::filesystem::exists(jpeg), "JPEG output should exist");
  require(std::filesystem::file_size(jpeg) > 0, "JPEG output should not be empty");

  std::filesystem::remove(tiff);
  std::filesystem::remove(jpeg);
}

void reads_linear_rgb_tiff_round_trip() {
  const astrocfa::RgbImage expected = tiny_image();
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "astrocfa-linear-rgb-roundtrip.tif";
  astrocfa::write_rgb_tiff16(expected, path.string());
  const astrocfa::RgbImage actual = astrocfa::read_linear_rgb_tiff(path.string());

  require(actual.width() == expected.width() && actual.height() == expected.height(),
          "TIFF round-trip dimensions");
  const astrocfa::RgbPixel pixel = actual.pixel(1, 0);
  require(std::abs(pixel.r - 0.3F) < 2.0F / 65535.0F,
          "TIFF reader should preserve linear red");
  require(std::abs(pixel.g - 0.4F) < 2.0F / 65535.0F,
          "TIFF reader should preserve linear green");
  require(std::abs(pixel.b - 0.5F) < 2.0F / 65535.0F,
          "TIFF reader should preserve linear blue");
  std::filesystem::remove(path);
}

void writes_loadable_linear_cfa_dng() {
  const astrocfa::BayerPattern pattern;
  astrocfa::CfaFrame expected(128, 96, pattern);
  for(std::size_t y = 0; y < expected.height(); ++y) {
    for(std::size_t x = 0; x < expected.width(); ++x) {
      expected.set_sample(x, y, static_cast<float>(x + y) / 256.0F);
    }
  }

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "astrocfa-linear-cfa-roundtrip.dng";
  astrocfa::write_linear_cfa_dng16(expected, path.string());
  const astrocfa::LinearRawFrame actual = astrocfa::load_linear_cfa_file(path.string());

  require(actual.cfa.width() == expected.width() &&
              actual.cfa.height() == expected.height(),
          "DNG round-trip dimensions");
  require(actual.cfa.pattern().at(0, 0) == astrocfa::CfaColor::red &&
              actual.cfa.pattern().at(1, 1) == astrocfa::CfaColor::blue,
          "DNG round-trip Bayer phase");
  require(std::abs(actual.cfa.sample(7, 5) - expected.sample(7, 5)) <
              2.0F / 65535.0F,
          "DNG round-trip normalized sample");
  std::filesystem::remove(path);
}

} // namespace

int main() {
  try {
    writes_tiff_and_jpeg();
    reads_linear_rgb_tiff_round_trip();
    writes_loadable_linear_cfa_dng();
  } catch(const std::exception &error) {
    std::cerr << "image_writer_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
