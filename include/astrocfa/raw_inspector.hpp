#pragma once

#include <array>
#include <cstdint>
#include <ostream>
#include <string>

namespace astrocfa {

struct PhaseStats {
  std::uint64_t samples = 0;
  std::uint64_t clipped = 0;
  std::uint16_t minimum = 0;
  std::uint16_t maximum = 0;
  double mean = 0.0;
};

struct RawInspection {
  std::string path;
  std::string make;
  std::string model;
  std::string normalized_make;
  std::string normalized_model;
  std::string cfa_description;

  int colors = 0;
  unsigned filters = 0;

  int raw_width = 0;
  int raw_height = 0;
  int active_width = 0;
  int active_height = 0;
  int top_margin = 0;
  int left_margin = 0;

  unsigned black_level = 0;
  std::array<unsigned, 4> cblack = {0, 0, 0, 0};
  unsigned white_level = 0;
  unsigned data_maximum = 0;
  unsigned raw_bits_per_sample = 0;

  float iso = 0.0F;
  float exposure_seconds = 0.0F;
  float aperture = 0.0F;
  float focal_length = 0.0F;

  bool has_bayer_cfa = false;
  bool has_xtrans_cfa = false;
  std::array<std::array<int, 2>, 2> bayer_phase = {{{0, 0}, {0, 0}}};
  std::array<PhaseStats, 4> phase_stats;

  std::uint64_t clipped_samples = 0;
  std::uint64_t clipped_components = 0;
  std::uint64_t largest_clipped_component = 0;

  double green_split_percent = 0.0;
};

RawInspection inspect_raw_file(const std::string &path);
void write_inspection_report(const RawInspection &inspection, std::ostream &out);

} // namespace astrocfa
