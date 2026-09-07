#include "astrocfa/frequency_cfa.hpp"

#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void uniform_signal_has_no_carrier_risk() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, astrocfa::CfaSample{.value = 0.5F, .valid = true, .clipped = false});
    }
  }

  const astrocfa::FrequencyCfaDiagnostics diagnostics =
      astrocfa::analyze_frequency_cfa(cfa, astrocfa::FrequencyCfaOptions{.tile_size = 4});
  require(diagnostics.total_tiles == 4, "Uniform frame should produce four tiles");
  require(diagnostics.max_alias_risk == 0.0, "Uniform signal should not carry Bayer modulation");
}

void checkerboard_signal_has_high_carrier_risk() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const float value = ((x + y) % 2U) == 0U ? 1.0F : 0.0F;
      cfa.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true, .clipped = false});
    }
  }

  const astrocfa::FrequencyCfaDiagnostics diagnostics =
      astrocfa::analyze_frequency_cfa(cfa, astrocfa::FrequencyCfaOptions{.tile_size = 4});
  require(diagnostics.high_risk_tiles == 4,
          "Checkerboard signal should trigger high CFA carrier risk");
  require(diagnostics.max_alias_risk > 0.9, "Checkerboard risk should be near one");
}

void risk_map_exposes_pixel_risk_from_tiles() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const float value = ((x + y) % 2U) == 0U ? 1.0F : 0.0F;
      cfa.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true, .clipped = false});
    }
  }

  const astrocfa::FrequencyRiskMap risk =
      astrocfa::build_frequency_risk_map(cfa, astrocfa::FrequencyCfaOptions{.tile_size = 4});
  require(risk.tiles_x() == 2 && risk.tiles_y() == 2, "8x8 frame should produce 2x2 risk tiles");
  require(risk.pixel_risk(6, 6) > 0.9, "Pixel risk should follow tile risk");
}

} // namespace

int main() {
  try {
    uniform_signal_has_no_carrier_risk();
    checkerboard_signal_has_high_carrier_risk();
    risk_map_exposes_pixel_risk_from_tiles();
  } catch(const std::exception &error) {
    std::cerr << "frequency_cfa_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
