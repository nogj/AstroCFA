#include "astrocfa/frequency_cfa.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

double carrier_sign(std::size_t x, std::size_t y, int carrier) {
  switch(carrier) {
  case 0:
    return (x % 2U) == 0U ? 1.0 : -1.0;
  case 1:
    return (y % 2U) == 0U ? 1.0 : -1.0;
  default:
    return ((x + y) % 2U) == 0U ? 1.0 : -1.0;
  }
}

} // namespace

namespace astrocfa {

FrequencyRiskMap::FrequencyRiskMap(std::size_t tiles_x, std::size_t tiles_y,
                                   std::size_t tile_size)
    : tiles_x_(tiles_x), tiles_y_(tiles_y), tile_size_(tile_size),
      risk_(tiles_x * tiles_y, 0.0) {
  if(tiles_x == 0 || tiles_y == 0 || tile_size == 0) {
    throw std::invalid_argument("Frequency risk map dimensions must be non-empty");
  }
}

double FrequencyRiskMap::tile_risk(std::size_t tile_x, std::size_t tile_y) const {
  return risk_[offset(tile_x, tile_y)];
}

void FrequencyRiskMap::set_tile_risk(std::size_t tile_x, std::size_t tile_y,
                                     double risk) {
  risk_[offset(tile_x, tile_y)] = std::clamp(risk, 0.0, 1.0);
}

double FrequencyRiskMap::pixel_risk(std::size_t x, std::size_t y) const {
  const std::size_t tile_x = std::min(tiles_x_ - 1, x / tile_size_);
  const std::size_t tile_y = std::min(tiles_y_ - 1, y / tile_size_);
  return tile_risk(tile_x, tile_y);
}

std::size_t FrequencyRiskMap::offset(std::size_t tile_x, std::size_t tile_y) const {
  if(tile_x >= tiles_x_ || tile_y >= tiles_y_) {
    throw std::out_of_range("Frequency risk tile coordinate out of range");
  }
  return tile_y * tiles_x_ + tile_x;
}

namespace {

double compute_tile_alias_risk(const CfaFrame &cfa, std::size_t x0, std::size_t y0,
                               std::size_t x1, std::size_t y1) {
  double sum = 0.0;
  std::size_t samples = 0;
  for(std::size_t y = y0; y < y1; ++y) {
    for(std::size_t x = x0; x < x1; ++x) {
      const CfaSample sample = cfa.sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }
      sum += sample.value;
      samples += 1;
    }
  }

  if(samples < 4) {
    return 0.0;
  }

  const double mean = sum / static_cast<double>(samples);
  double ac_energy = 0.0;
  double carrier_projection[3] = {0.0, 0.0, 0.0};

  for(std::size_t y = y0; y < y1; ++y) {
    for(std::size_t x = x0; x < x1; ++x) {
      const CfaSample sample = cfa.sample_info(x, y);
      if(!sample.valid || sample.clipped) {
        continue;
      }

      const double centered = sample.value - mean;
      ac_energy += centered * centered;
      for(int carrier = 0; carrier < 3; ++carrier) {
        carrier_projection[carrier] += centered * carrier_sign(x, y, carrier);
      }
    }
  }

  ac_energy /= static_cast<double>(samples);
  double carrier_energy = 0.0;
  for(double projection : carrier_projection) {
    const double normalized_projection = projection / static_cast<double>(samples);
    carrier_energy += normalized_projection * normalized_projection;
  }

  return ac_energy > 1.0e-12 ? std::min(1.0, carrier_energy / ac_energy) : 0.0;
}

} // namespace

FrequencyCfaDiagnostics analyze_frequency_cfa(const CfaFrame &cfa,
                                              FrequencyCfaOptions options) {
  if(options.tile_size == 0) {
    throw std::invalid_argument("Frequency CFA tile size must be greater than zero");
  }

  FrequencyCfaDiagnostics diagnostics;
  diagnostics.tile_size = options.tile_size;

  for(std::size_t y0 = 0; y0 < cfa.height(); y0 += options.tile_size) {
    for(std::size_t x0 = 0; x0 < cfa.width(); x0 += options.tile_size) {
      const std::size_t x1 = std::min(cfa.width(), x0 + options.tile_size);
      const std::size_t y1 = std::min(cfa.height(), y0 + options.tile_size);

      double sum = 0.0;
      std::size_t samples = 0;
      double ac_energy = 0.0;
      double carrier_projection[3] = {0.0, 0.0, 0.0};
      for(std::size_t y = y0; y < y1; ++y) {
        for(std::size_t x = x0; x < x1; ++x) {
          const CfaSample sample = cfa.sample_info(x, y);
          if(!sample.valid || sample.clipped) {
            continue;
          }
          sum += sample.value;
          samples += 1;
        }
      }

      if(samples < 4) {
        continue;
      }

      const double mean = sum / static_cast<double>(samples);

      for(std::size_t y = y0; y < y1; ++y) {
        for(std::size_t x = x0; x < x1; ++x) {
          const CfaSample sample = cfa.sample_info(x, y);
          if(!sample.valid || sample.clipped) {
            continue;
          }
          const double centered = sample.value - mean;
          ac_energy += centered * centered;
          for(int carrier = 0; carrier < 3; ++carrier) {
            carrier_projection[carrier] += centered * carrier_sign(x, y, carrier);
          }
        }
      }

      ac_energy /= static_cast<double>(samples);
      double carrier_energy = 0.0;
      for(double projection : carrier_projection) {
        const double normalized_projection = projection / static_cast<double>(samples);
        carrier_energy += normalized_projection * normalized_projection;
      }

      const double alias_risk =
          ac_energy > 1.0e-12 ? std::min(1.0, carrier_energy / ac_energy) : 0.0;

      diagnostics.total_tiles += 1;
      diagnostics.high_risk_tiles +=
          alias_risk >= options.high_risk_threshold ? 1U : 0U;
      diagnostics.mean_ac_energy += ac_energy;
      diagnostics.mean_carrier_energy += carrier_energy;
      diagnostics.mean_alias_risk += alias_risk;
      diagnostics.max_alias_risk = std::max(diagnostics.max_alias_risk, alias_risk);
    }
  }

  if(diagnostics.total_tiles > 0) {
    const double tiles = static_cast<double>(diagnostics.total_tiles);
    diagnostics.mean_ac_energy /= tiles;
    diagnostics.mean_carrier_energy /= tiles;
    diagnostics.mean_alias_risk /= tiles;
  }

  return diagnostics;
}

FrequencyRiskMap build_frequency_risk_map(const CfaFrame &cfa, FrequencyCfaOptions options) {
  if(options.tile_size == 0) {
    throw std::invalid_argument("Frequency CFA tile size must be greater than zero");
  }

  const std::size_t tiles_x = (cfa.width() + options.tile_size - 1U) / options.tile_size;
  const std::size_t tiles_y = (cfa.height() + options.tile_size - 1U) / options.tile_size;
  FrequencyRiskMap risk_map(tiles_x, tiles_y, options.tile_size);

  for(std::size_t tile_y = 0; tile_y < tiles_y; ++tile_y) {
    for(std::size_t tile_x = 0; tile_x < tiles_x; ++tile_x) {
      const std::size_t x0 = tile_x * options.tile_size;
      const std::size_t y0 = tile_y * options.tile_size;
      const std::size_t x1 = std::min(cfa.width(), x0 + options.tile_size);
      const std::size_t y1 = std::min(cfa.height(), y0 + options.tile_size);
      risk_map.set_tile_risk(tile_x, tile_y, compute_tile_alias_risk(cfa, x0, y0, x1, y1));
    }
  }

  return risk_map;
}

} // namespace astrocfa
