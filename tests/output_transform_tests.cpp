#include "astrocfa/output_transform.hpp"

#include <cassert>

int main() {
  astrocfa::RgbImage image(4, 1);
  image.set_pixel(0, 0, astrocfa::RgbPixel{.r = 0.0F, .g = 0.0F, .b = 0.0F});
  image.set_pixel(1, 0, astrocfa::RgbPixel{.r = 0.01F, .g = 0.01F, .b = 0.01F});
  image.set_pixel(2, 0, astrocfa::RgbPixel{.r = 0.10F, .g = 0.10F, .b = 0.10F});
  image.set_pixel(3, 0, astrocfa::RgbPixel{.r = 1.0F, .g = 1.0F, .b = 1.0F});

  const astrocfa::RgbImage preview = astrocfa::make_astro_preview(image);
  assert(preview.width() == image.width());
  assert(preview.height() == image.height());
  assert(preview.pixel(0, 0).g <= preview.pixel(1, 0).g);
  assert(preview.pixel(1, 0).g <= preview.pixel(2, 0).g);
  assert(preview.pixel(2, 0).g <= preview.pixel(3, 0).g);
  assert(preview.pixel(3, 0).g <= 1.0F);

  return 0;
}
