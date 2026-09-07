#include "astrocfa/raw_inspector.hpp"

#include "astrocfa/connected_components.hpp"

#include <libraw/libraw.h>

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

std::string c_string(const char *value) {
  return value == nullptr ? std::string{} : std::string(value);
}

std::string color_name(const astrocfa::RawInspection &inspection, int color) {
  if(color >= 0 && color < 4 && inspection.cfa_description.size() > static_cast<std::size_t>(color)) {
    return std::string(1, inspection.cfa_description[static_cast<std::size_t>(color)]);
  }
  return "C" + std::to_string(color);
}

std::string phase_name(const astrocfa::RawInspection &inspection, int color) {
  return "C" + std::to_string(color) + "/" + color_name(inspection, color);
}

} // namespace

namespace astrocfa {

RawInspection inspect_raw_file(const std::string &path) {
  LibRaw raw;
  const int open_result = raw.open_file(path.c_str());
  if(open_result != LIBRAW_SUCCESS) {
    throw std::runtime_error("Cannot open RAW file: " + std::string(libraw_strerror(open_result)));
  }

  const int unpack_result = raw.unpack();
  if(unpack_result != LIBRAW_SUCCESS) {
    throw std::runtime_error("Cannot unpack RAW file: " + std::string(libraw_strerror(unpack_result)));
  }

  const auto &idata = raw.imgdata.idata;
  const auto &sizes = raw.imgdata.sizes;
  const auto &color = raw.imgdata.color;
  const auto &other = raw.imgdata.other;

  RawInspection out;
  out.path = path;
  out.make = c_string(idata.make);
  out.model = c_string(idata.model);
  out.normalized_make = c_string(idata.normalized_make);
  out.normalized_model = c_string(idata.normalized_model);
  out.cfa_description = c_string(idata.cdesc);
  out.colors = idata.colors;
  out.filters = idata.filters;
  out.raw_width = sizes.raw_width;
  out.raw_height = sizes.raw_height;
  out.active_width = sizes.width;
  out.active_height = sizes.height;
  out.top_margin = sizes.top_margin;
  out.left_margin = sizes.left_margin;
  out.black_level = color.black;
  out.white_level = color.maximum;
  out.data_maximum = color.data_maximum;
  out.raw_bits_per_sample = color.raw_bps;
  out.iso = other.iso_speed;
  out.exposure_seconds = other.shutter;
  out.aperture = other.aperture;
  out.focal_length = other.focal_len;

  for(std::size_t i = 0; i < out.cblack.size(); ++i) {
    out.cblack[i] = color.cblack[i];
  }

  out.has_bayer_cfa = idata.filters != 0 && idata.filters >= 1000 && idata.colors >= 3;
  out.has_xtrans_cfa = idata.filters == 9;

  if(!raw.imgdata.rawdata.raw_image) {
    return out;
  }

  for(auto &stats : out.phase_stats) {
    stats.minimum = std::numeric_limits<std::uint16_t>::max();
  }

  const auto white = static_cast<std::uint16_t>(
      std::min<unsigned>(out.white_level == 0 ? std::numeric_limits<std::uint16_t>::max()
                                              : out.white_level,
                         std::numeric_limits<std::uint16_t>::max()));
  std::vector<std::uint8_t> clipped_mask(static_cast<std::size_t>(out.active_width) *
                                             static_cast<std::size_t>(out.active_height),
                                         0);

  for(int y = 0; y < out.active_height; ++y) {
    const int raw_y = y + out.top_margin;
    if(raw_y < 0 || raw_y >= out.raw_height) {
      continue;
    }

    for(int x = 0; x < out.active_width; ++x) {
      const int raw_x = x + out.left_margin;
      if(raw_x < 0 || raw_x >= out.raw_width) {
        continue;
      }

      const int color_index = raw.COLOR(raw_y, raw_x);
      if(color_index < 0 || color_index >= 4) {
        continue;
      }

      if(y < 2 && x < 2) {
        out.bayer_phase[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)] =
            color_index;
      }

      const auto value = raw.imgdata.rawdata
                             .raw_image[static_cast<std::size_t>(raw_y) *
                                            static_cast<std::size_t>(out.raw_width) +
                                        static_cast<std::size_t>(raw_x)];
      auto &stats = out.phase_stats[static_cast<std::size_t>(color_index)];
      stats.samples += 1;
      const bool clipped = value >= white;
      stats.clipped += clipped ? 1 : 0;
      if(clipped) {
        clipped_mask[static_cast<std::size_t>(y) * static_cast<std::size_t>(out.active_width) +
                     static_cast<std::size_t>(x)] = 1;
      }
      stats.minimum = std::min(stats.minimum, value);
      stats.maximum = std::max(stats.maximum, value);
      stats.mean += static_cast<double>(value);
    }
  }

  for(auto &stats : out.phase_stats) {
    if(stats.samples > 0) {
      stats.mean /= static_cast<double>(stats.samples);
    } else {
      stats.minimum = 0;
    }
  }

  const auto &green1 = out.phase_stats[1];
  const auto &green2 = out.phase_stats[3];
  if(green1.samples > 0 && green2.samples > 0) {
    const double green_average = 0.5 * (green1.mean + green2.mean);
    if(green_average > 0.0) {
      out.green_split_percent = 100.0 * (green1.mean - green2.mean) / green_average;
    }
  }

  const BinaryComponentStats clipped_components = analyze_binary_components(
      clipped_mask, static_cast<std::size_t>(out.active_width),
      static_cast<std::size_t>(out.active_height));
  out.clipped_samples = clipped_components.marked_samples;
  out.clipped_components = clipped_components.components;
  out.largest_clipped_component = clipped_components.largest_area;

  return out;
}

void write_inspection_report(const RawInspection &inspection, std::ostream &out) {
  out << "AstroCFA RAW inspection\n"
      << "  file: " << inspection.path << "\n"
      << "  camera: " << inspection.normalized_make << " " << inspection.normalized_model
      << "\n"
      << "  original camera: " << inspection.make << " " << inspection.model << "\n"
      << "  raw size: " << inspection.raw_width << " x " << inspection.raw_height << "\n"
      << "  active area: " << inspection.active_width << " x " << inspection.active_height
      << " @ +" << inspection.left_margin << ",+" << inspection.top_margin << "\n"
      << "  cfa: " << (inspection.has_bayer_cfa ? "Bayer" : inspection.has_xtrans_cfa ? "X-Trans" : "unknown")
      << "  colors=" << inspection.colors << "  cdesc=" << inspection.cfa_description
      << "  filters=0x" << std::hex << inspection.filters << std::dec << "\n"
      << "  black: " << inspection.black_level << "  cblack: [" << inspection.cblack[0]
      << ", " << inspection.cblack[1] << ", " << inspection.cblack[2] << ", "
      << inspection.cblack[3] << "]\n"
      << "  white: " << inspection.white_level << "  data maximum: "
      << inspection.data_maximum << "  raw bits: " << inspection.raw_bits_per_sample << "\n"
      << "  exposure: ISO " << inspection.iso << ", " << inspection.exposure_seconds
      << " s, f/" << inspection.aperture << ", " << inspection.focal_length << " mm\n";

  if(inspection.has_bayer_cfa) {
    out << "  bayer phase 2x2: "
        << color_name(inspection, inspection.bayer_phase[0][0])
        << color_name(inspection, inspection.bayer_phase[0][1]) << " / "
        << color_name(inspection, inspection.bayer_phase[1][0])
        << color_name(inspection, inspection.bayer_phase[1][1]) << "\n";
  }

  out << "\nPhase diagnostics over active area:\n"
      << "  phase  samples       min    mean      max    clipped\n";

  out << std::fixed << std::setprecision(2);
  for(std::size_t i = 0; i < inspection.phase_stats.size(); ++i) {
    const auto &stats = inspection.phase_stats[i];
    if(stats.samples == 0) {
      continue;
    }

    const double clipped_percent =
        100.0 * static_cast<double>(stats.clipped) / static_cast<double>(stats.samples);
    out << "  " << std::left << std::setw(5)
        << phase_name(inspection, static_cast<int>(i)) << std::right
        << std::setw(10) << stats.samples << "  " << std::setw(6) << stats.minimum
        << "  " << std::setw(8) << stats.mean << "  " << std::setw(6) << stats.maximum
        << "  " << stats.clipped << " (" << clipped_percent << "%)\n";
  }

  out << "\nSaturation topology:\n"
      << "  clipped samples: " << inspection.clipped_samples << "\n"
      << "  clipped components: " << inspection.clipped_components << "\n"
      << "  largest clipped component: " << inspection.largest_clipped_component
      << " samples\n";

  if(inspection.phase_stats[1].samples > 0 && inspection.phase_stats[3].samples > 0) {
    out << "\nCFA balance diagnostics:\n"
        << "  G1/G2 mean split: " << inspection.green_split_percent << "%\n";
  }

  out << "\nAstro notes:\n";
  out << "  - clipped phase counts identify saturated stars before demosaicing.\n";
  out << "  - per-phase means are useful for flat-field and CFA balance sanity checks.\n";
  out << "  - remosaicing residual maps will build on these same active-area coordinates.\n";
}

} // namespace astrocfa
