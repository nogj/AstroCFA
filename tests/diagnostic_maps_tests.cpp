#include "astrocfa/diagnostic_maps.hpp"
#include "astrocfa/demosaic.hpp"

#include <cassert>

int main() {
  astrocfa::CfaFrame cfa(8, 8, astrocfa::BayerPattern{});
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, static_cast<float>((x + y) % 4U) / 4.0F);
    }
  }
  cfa.set_clipped(0, 0, true);

  const astrocfa::RgbImage alias =
      astrocfa::make_frequency_alias_risk_map(cfa, astrocfa::FrequencyCfaOptions{
                                                       .tile_size = 4,
                                                       .high_risk_threshold = 0.35,
                                                   });
  assert(alias.width() == cfa.width());
  assert(alias.height() == cfa.height());

  const astrocfa::RgbImage reconstruction = astrocfa::demosaic_bilinear_baseline(cfa);
  const astrocfa::RgbImage residual =
      astrocfa::make_remosaic_residual_map(cfa, reconstruction);
  assert(residual.width() == cfa.width());
  assert(residual.height() == cfa.height());
  assert(residual.pixel(0, 0).r == 1.0F);
  assert(residual.pixel(0, 0).g == 0.0F);

  astrocfa::DefectMap defects(8, 8);
  defects.add(2, 2, astrocfa::SensorDefect::hot);
  defects.add(3, 3, astrocfa::SensorDefect::dead);
  defects.add(4, 4, astrocfa::SensorDefect::invalid_master);
  const astrocfa::RgbImage defect_image = astrocfa::make_sensor_defect_map(defects);
  assert(defect_image.pixel(2, 2).r == 1.0F);
  assert(defect_image.pixel(3, 3).b == 1.0F);
  assert(defect_image.pixel(4, 4).g == 1.0F);

  return 0;
}
