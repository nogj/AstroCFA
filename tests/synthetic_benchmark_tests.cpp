#include "astrocfa/demosaic.hpp"
#include "astrocfa/reconstruction_metrics.hpp"
#include "astrocfa/synthetic_astro_scene.hpp"

#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void synthetic_scene_has_truth_and_cfa() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 48,
          .height = 32,
          .seed = 7,
          .add_noise = false,
          .add_hot_pixels = false,
      });
  require(scene.truth.width() == 48, "Synthetic truth width");
  require(scene.cfa.height() == 32, "Synthetic CFA height");
  require(!scene.stars.empty(), "Synthetic scene should include stars");
}

void benchmark_metrics_are_finite_and_accountable() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 64,
          .height = 48,
          .seed = 11,
          .add_noise = false,
          .add_hot_pixels = false,
      });
  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_inverse_refine(scene.cfa, astrocfa::NoiseModel{},
                                           astrocfa::InverseRefinementOptions{
                                               .frequency = {.tile_size = 8},
                                               .iterations = 2,
                                           });
  const astrocfa::ReconstructionMetrics metrics =
      astrocfa::measure_reconstruction(scene.truth, result.image, scene.cfa, scene.stars);
  require(metrics.samples == scene.truth.width() * scene.truth.height() * 3U,
          "Metric should count RGB samples");
  require(metrics.rgb_mae >= 0.0, "RGB MAE should be finite and non-negative");
  require(metrics.star_samples > 0, "Star metrics should cover star pixels");
  require(metrics.cfa_residual_mae == 0.0,
          "Measurement-preserving reconstruction should have zero CFA MAE");
}

} // namespace

int main() {
  try {
    synthetic_scene_has_truth_and_cfa();
    benchmark_metrics_are_finite_and_accountable();
  } catch(const std::exception &error) {
    std::cerr << "synthetic_benchmark_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
