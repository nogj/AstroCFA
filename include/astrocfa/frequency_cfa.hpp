#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>
#include <vector>

namespace astrocfa {

struct FrequencyCfaOptions {
  std::size_t tile_size = 64;
  double high_risk_threshold = 0.35;
};

struct FrequencyCfaDiagnostics {
  std::size_t tile_size = 0;
  std::size_t total_tiles = 0;
  std::size_t high_risk_tiles = 0;
  double mean_ac_energy = 0.0;
  double mean_carrier_energy = 0.0;
  double mean_alias_risk = 0.0;
  double max_alias_risk = 0.0;
};

class FrequencyRiskMap {
public:
  FrequencyRiskMap(std::size_t tiles_x, std::size_t tiles_y, std::size_t tile_size);

  [[nodiscard]] std::size_t tiles_x() const { return tiles_x_; }
  [[nodiscard]] std::size_t tiles_y() const { return tiles_y_; }
  [[nodiscard]] std::size_t tile_size() const { return tile_size_; }

  [[nodiscard]] double tile_risk(std::size_t tile_x, std::size_t tile_y) const;
  void set_tile_risk(std::size_t tile_x, std::size_t tile_y, double risk);
  [[nodiscard]] double pixel_risk(std::size_t x, std::size_t y) const;

private:
  [[nodiscard]] std::size_t offset(std::size_t tile_x, std::size_t tile_y) const;

  std::size_t tiles_x_ = 0;
  std::size_t tiles_y_ = 0;
  std::size_t tile_size_ = 1;
  std::vector<double> risk_;
};

[[nodiscard]] FrequencyCfaDiagnostics
analyze_frequency_cfa(const CfaFrame &cfa, FrequencyCfaOptions options = {});
[[nodiscard]] FrequencyRiskMap build_frequency_risk_map(const CfaFrame &cfa,
                                                        FrequencyCfaOptions options = {});

} // namespace astrocfa
