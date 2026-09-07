#include "astrocfa/demosaic.hpp"
#include "astrocfa/multiframe_benchmark.hpp"
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
  require(metrics.measured_stars == scene.stars.size(),
          "Photometry and shape metrics should cover every synthetic star");
  require(std::isfinite(metrics.star_flux_relative_error),
          "Star photometry error should be finite");
  require(std::isfinite(metrics.star_fwhm_relative_error),
          "Star FWHM error should be finite");
  require(std::isfinite(metrics.star_elongation_error),
          "Star elongation error should be finite");
  require(metrics.cfa_residual_mae == 0.0,
          "Measurement-preserving reconstruction should have zero CFA MAE");
}

void star_chroma_guard_reduces_synthetic_star_false_color() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 96,
          .height = 72,
          .seed = 5,
          .add_noise = true,
          .add_hot_pixels = true,
      });
  const astrocfa::InverseRefinementOptions no_guard{
      .frequency = {.tile_size = 16},
      .iterations = 3,
      .star_chroma_guard = 0.0,
  };
  const astrocfa::InverseRefinementOptions guarded{
      .frequency = {.tile_size = 16},
      .iterations = 3,
      .star_chroma_guard = 0.35,
  };
  const astrocfa::DemosaicResult unguarded =
      astrocfa::reconstruct_inverse_refine(scene.cfa, astrocfa::NoiseModel{}, no_guard);
  const astrocfa::DemosaicResult guarded_result =
      astrocfa::reconstruct_inverse_refine(scene.cfa, astrocfa::NoiseModel{}, guarded);
  const astrocfa::ReconstructionMetrics unguarded_metrics =
      astrocfa::measure_reconstruction(scene.truth, unguarded.image, scene.cfa,
                                       scene.stars);
  const astrocfa::ReconstructionMetrics guarded_metrics =
      astrocfa::measure_reconstruction(scene.truth, guarded_result.image, scene.cfa,
                                       scene.stars);

  require(guarded_result.residual.maximum_absolute == 0.0F,
          "Star chroma guard must preserve measured CFA samples");
  require(guarded_metrics.star_false_color <= unguarded_metrics.star_false_color,
          "Star chroma guard should reduce synthetic star false color");
}

void multiframe_benchmark_is_reproducible_and_auditable() {
  const astrocfa::MultiframeBenchmarkResult benchmark =
      astrocfa::run_multiframe_benchmark(astrocfa::MultiframeBenchmarkOptions{
          .width = 48,
          .height = 32,
          .frames = 4,
          .iterations = 2,
          .transients_per_frame = 2,
          .seed = 17,
          .add_noise = true,
          .vary_seeing = false,
      });
  require(benchmark.methods.size() == 4, "Benchmark should compare four pipelines");
  require(benchmark.injected_transients == 8,
          "Benchmark should report every injected transient");
  require(benchmark.methods.back().name == "joint-robust",
          "Robust joint solver should be the final benchmark method");
  require(benchmark.methods.back().solver_stats.robust_outliers > 0,
          "Robust solver should expose inconsistent measurements");
  const auto &individual = benchmark.methods.front().metrics;
  const auto &robust = benchmark.methods.back().metrics;
  require(robust.rgb_rmse < individual.rgb_rmse,
          "Joint reconstruction should beat individual demosaic RGB error");
  require(robust.star_luma_rmse < individual.star_luma_rmse,
          "Joint reconstruction should improve stellar luminance error");
  require(robust.star_fwhm_relative_error < individual.star_fwhm_relative_error,
          "Joint reconstruction should improve stellar FWHM recovery");
  for(const auto &method : benchmark.methods) {
    require(std::isfinite(method.metrics.rgb_rmse),
            "Every multiframe benchmark metric should be finite");
    require(method.metrics.measured_stars == benchmark.stars.size(),
            "Every multiframe method should measure all stars");
  }
}

} // namespace

int main() {
  try {
    synthetic_scene_has_truth_and_cfa();
    benchmark_metrics_are_finite_and_accountable();
    star_chroma_guard_reduces_synthetic_star_false_color();
    multiframe_benchmark_is_reproducible_and_auditable();
  } catch(const std::exception &error) {
    std::cerr << "synthetic_benchmark_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
