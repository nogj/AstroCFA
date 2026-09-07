#include "astrocfa/image_writer.hpp"

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

} // namespace

int main() {
  try {
    writes_tiff_and_jpeg();
  } catch(const std::exception &error) {
    std::cerr << "image_writer_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}

