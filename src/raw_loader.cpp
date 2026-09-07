#include "astrocfa/raw_loader.hpp"

#include <libraw/libraw.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

astrocfa::BayerPattern active_bayer_pattern(LibRaw &raw) {
  const auto &sizes = raw.imgdata.sizes;
  astrocfa::BayerPattern pattern;

  for(std::size_t y = 0; y < 2; ++y) {
    for(std::size_t x = 0; x < 2; ++x) {
      const int color = raw.COLOR(static_cast<int>(sizes.top_margin + y),
                                  static_cast<int>(sizes.left_margin + x));
      pattern.phase[y * 2U + x] = static_cast<astrocfa::CfaColor>(std::clamp(color, 0, 3));
    }
  }

  return pattern;
}

float normalized_linear_value(unsigned raw_value, unsigned black, unsigned white) {
  if(white <= black) {
    return 0.0F;
  }

  const float numerator =
      static_cast<float>(std::max<int>(0, static_cast<int>(raw_value) - static_cast<int>(black)));
  const float denominator = static_cast<float>(white - black);
  return numerator / denominator;
}

} // namespace

namespace astrocfa {

LinearRawFrame load_linear_cfa_file(const std::string &path) {
  LibRaw raw;
  const int open_result = raw.open_file(path.c_str());
  if(open_result != LIBRAW_SUCCESS) {
    throw std::runtime_error("Cannot open RAW file: " + std::string(libraw_strerror(open_result)));
  }

  const int unpack_result = raw.unpack();
  if(unpack_result != LIBRAW_SUCCESS) {
    throw std::runtime_error("Cannot unpack RAW file: " + std::string(libraw_strerror(unpack_result)));
  }

  if(!raw.imgdata.rawdata.raw_image) {
    throw std::runtime_error("RAW file does not expose an unpacked CFA plane");
  }

  const auto &idata = raw.imgdata.idata;
  if(!(idata.filters != 0 && idata.filters >= 1000 && idata.colors >= 3)) {
    throw std::runtime_error("Only Bayer CFA RAW files are supported by the linear loader");
  }

  LinearRawFrame frame{
      .inspection = inspect_raw_file(path),
      .cfa = CfaFrame(static_cast<std::size_t>(raw.imgdata.sizes.width),
                      static_cast<std::size_t>(raw.imgdata.sizes.height),
                      active_bayer_pattern(raw)),
  };

  const auto &sizes = raw.imgdata.sizes;
  const auto &color = raw.imgdata.color;
  const unsigned fallback_white =
      color.maximum > 0 ? color.maximum : std::numeric_limits<std::uint16_t>::max();

  for(std::size_t y = 0; y < frame.cfa.height(); ++y) {
    const auto raw_y = static_cast<std::size_t>(sizes.top_margin) + y;
    for(std::size_t x = 0; x < frame.cfa.width(); ++x) {
      const auto raw_x = static_cast<std::size_t>(sizes.left_margin) + x;
      const int phase_index = raw.COLOR(static_cast<int>(raw_y), static_cast<int>(raw_x));
      const unsigned black =
          color.black + (phase_index >= 0 && phase_index < 4 ? color.cblack[phase_index] : 0U);
      const unsigned raw_value =
          raw.imgdata.rawdata.raw_image[raw_y * static_cast<std::size_t>(sizes.raw_width) + raw_x];

      CfaSample sample;
      sample.clipped = raw_value >= fallback_white;
      sample.valid = true;
      sample.value = normalized_linear_value(raw_value, black, fallback_white);
      frame.cfa.set_sample(x, y, sample);
    }
  }

  return frame;
}

} // namespace astrocfa

