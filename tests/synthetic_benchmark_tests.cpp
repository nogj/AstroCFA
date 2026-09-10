#include "astrocfa/demosaic.hpp"
#include "astrocfa/multiframe_benchmark.hpp"
#include "astrocfa/morphological_reconstruction.hpp"
#include "astrocfa/reconstruction_metrics.hpp"
#include "astrocfa/synthetic_astro_scene.hpp"

#include <algorithm>
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
  require(std::isfinite(metrics.star_flux_relative_bias),
          "Star photometry bias should be finite");
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

void morphological_poc_is_sensor_domain_and_finite() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 96,
          .height = 72,
          .seed = 7,
          .add_noise = true,
          .add_hot_pixels = true,
      });
  const astrocfa::DemosaicResult result = astrocfa::reconstruct_morphological_cfa(
      scene.cfa,
      astrocfa::NoiseModel{.read_noise = 0.0025, .shot_noise_scale = 0.0018});
  const astrocfa::ReconstructionMetrics metrics = astrocfa::measure_reconstruction(
      scene.truth, result.image, scene.cfa, scene.stars);
  require(std::isfinite(metrics.rgb_rmse),
          "Morphological reconstruction RGB error should be finite");
  require(std::isfinite(metrics.star_flux_relative_error),
          "Morphological reconstruction photometry should be finite");
  require(metrics.cfa_residual_mae > 0.0,
          "Sensor-domain likelihood must not hard-restore noisy CFA samples");
}

void shared_profile_improves_non_gaussian_psf_shape() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 128,
          .height = 96,
          .seed = 19,
          .add_noise = true,
          .add_hot_pixels = true,
          .common_star_sigma = 0.65,
          .star_flux_scale = 0.6,
          .redundant_star_phase = false,
          .star_count = 8,
          .moffat_beta = 2.5,
      });
  const astrocfa::NoiseModel noise{
      .read_noise = 0.0025,
      .shot_noise_scale = 0.0018,
  };
  const astrocfa::DemosaicResult independent =
      astrocfa::reconstruct_morphological_cfa(
          scene.cfa, noise,
          astrocfa::MorphologicalReconstructionOptions{
              .share_psf_across_sources = false,
          });
  const astrocfa::MorphologicalReconstructionResult shared =
      astrocfa::reconstruct_morphological_cfa_detailed(scene.cfa, noise);
  const astrocfa::ReconstructionMetrics independent_metrics =
      astrocfa::measure_reconstruction(scene.truth, independent.image, scene.cfa,
                                       scene.stars);
  const astrocfa::ReconstructionMetrics shared_metrics =
      astrocfa::measure_reconstruction(scene.truth, shared.reconstruction.image,
                                       scene.cfa, scene.stars);
  require(scene.stars.size() == 8,
          "Phase-diversity fixture should contain the requested stars");
  require(shared_metrics.star_false_color < independent_metrics.star_false_color,
          "Observable shared profile should reduce stellar false color");
  require(shared_metrics.star_fwhm_relative_error <
              independent_metrics.star_fwhm_relative_error,
          "Observable shared profile should improve non-Gaussian PSF shape");
  require(shared.stats.selected_model == astrocfa::MorphologicalModel::shared_profile,
          "CFA holdout should select the shared profile on a common Moffat PSF");
  require(shared.stats.profile_sources > 1 &&
              shared.stats.profile_information_gain > 0.0,
          "Profile selection should expose positive full-matrix information gain");
  require(shared.stats.profile_validation < shared.stats.independent_validation &&
              shared.stats.profile_validation < shared.stats.shared_validation,
          "Selected shared profile should win on unseen CFA samples");
}

void shared_epsf_recovers_anisotropic_psf() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 192,
          .height = 144,
          .seed = 7,
          .add_noise = true,
          .add_hot_pixels = true,
          .common_star_sigma = 0.65,
          .star_flux_scale = 0.5,
          .redundant_star_phase = false,
          .star_count = 16,
          .moffat_beta = 2.5,
          .psf_ellipticity = 0.35,
          .psf_angle = 0.55,
      });
  const astrocfa::NoiseModel noise{
      .read_noise = 0.0025,
      .shot_noise_scale = 0.0018,
  };
  const astrocfa::DemosaicResult without_epsf =
      astrocfa::reconstruct_morphological_cfa(
          scene.cfa, noise,
          astrocfa::MorphologicalReconstructionOptions{.enable_epsf = false});
  const astrocfa::MorphologicalReconstructionResult with_epsf =
      astrocfa::reconstruct_morphological_cfa_detailed(scene.cfa, noise);
  const astrocfa::ReconstructionMetrics baseline_metrics =
      astrocfa::measure_reconstruction(scene.truth, without_epsf.image, scene.cfa,
                                       scene.stars);
  const astrocfa::ReconstructionMetrics epsf_metrics =
      astrocfa::measure_reconstruction(scene.truth, with_epsf.reconstruction.image,
                                       scene.cfa, scene.stars);
  require(with_epsf.stats.selected_model == astrocfa::MorphologicalModel::shared_epsf,
          "CFA holdout should select a 2D ePSF for an anisotropic common PSF");
  require(with_epsf.stats.epsf_validation < with_epsf.stats.profile_validation,
          "The 2D ePSF should predict unseen anisotropic CFA samples better");
  require(with_epsf.stats.reconstructed_sources >
              with_epsf.stats.validated_sources,
          "The ePSF pass should recover sources rejected by the circular seed model");
  require(epsf_metrics.star_fwhm_relative_error <
              baseline_metrics.star_fwhm_relative_error,
          "The selected ePSF should improve anisotropic stellar FWHM");
  require(epsf_metrics.star_elongation_error <
              baseline_metrics.star_elongation_error,
          "The selected ePSF should improve anisotropic stellar shape");
}

void chromatic_epsf_reduces_lateral_color_error() {
  const astrocfa::SyntheticAstroScene scene =
      astrocfa::make_synthetic_astro_scene(astrocfa::SyntheticAstroSceneOptions{
          .width = 256,
          .height = 192,
          .seed = 7,
          .add_noise = true,
          .add_hot_pixels = true,
          .common_star_sigma = 0.65,
          .star_flux_scale = 0.6,
          .redundant_star_phase = false,
          .star_count = 20,
          .moffat_beta = 2.5,
          .psf_ellipticity = 0.20,
          .psf_angle = 0.55,
          .chromatic_psf_shift = 0.25,
          .chromatic_psf_scale = 0.10,
      });
  const astrocfa::NoiseModel noise{
      .read_noise = 0.0025,
      .shot_noise_scale = 0.0018,
  };
  const astrocfa::DemosaicResult achromatic =
      astrocfa::reconstruct_morphological_cfa(
          scene.cfa, noise,
          astrocfa::MorphologicalReconstructionOptions{
              .enable_chromatic_epsf = false,
          });
  const astrocfa::MorphologicalReconstructionResult chromatic =
      astrocfa::reconstruct_morphological_cfa_detailed(scene.cfa, noise);
  const astrocfa::ReconstructionMetrics achromatic_metrics =
      astrocfa::measure_reconstruction(scene.truth, achromatic.image, scene.cfa,
                                       scene.stars);
  const astrocfa::ReconstructionMetrics chromatic_metrics =
      astrocfa::measure_reconstruction(scene.truth,
                                       chromatic.reconstruction.image, scene.cfa,
                                       scene.stars);
  require(chromatic.stats.selected_model ==
              astrocfa::MorphologicalModel::shared_chromatic_epsf,
          "CFA holdout should select the low-rank chromatic ePSF");
  require(chromatic.stats.chromatic_epsf_validation <
              chromatic.stats.epsf_validation,
          "Chromatic ePSF should predict unseen CFA samples better");
  require(chromatic_metrics.star_false_color <
              achromatic_metrics.star_false_color,
          "Chromatic ePSF should reduce lateral stellar color error");
  require(chromatic_metrics.star_luma_rmse < achromatic_metrics.star_luma_rmse,
          "Chromatic ePSF should improve stellar luminance reconstruction");
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
  require(benchmark.methods.size() == 5, "Benchmark should compare five pipelines");
  require(benchmark.injected_transients == 8,
          "Benchmark should report every injected transient");
  require(benchmark.methods.back().name == "joint-robust-psf",
          "PSF-aware robust solver should be the final benchmark method");
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

void psf_aware_solver_recovers_variable_seeing_detail() {
  const astrocfa::MultiframeBenchmarkResult benchmark =
      astrocfa::run_multiframe_benchmark(astrocfa::MultiframeBenchmarkOptions{
          .width = 64,
          .height = 48,
          .frames = 4,
          .iterations = 5,
          .transients_per_frame = 0,
          .seed = 23,
          .add_noise = false,
          .vary_seeing = true,
      });
  const auto &individual = benchmark.methods.front().metrics;
  const auto &without_psf = benchmark.methods[3].metrics;
  const auto &with_psf = benchmark.methods.back().metrics;
  const auto auto_psf = std::find_if(
      benchmark.methods.begin(), benchmark.methods.end(),
      [](const astrocfa::MultiframeBenchmarkMethod &method) {
        return method.name == "joint-robust-auto-psf";
      });
  require(benchmark.methods[3].name == "joint-robust-no-psf",
          "Benchmark should retain the PSF ablation");
  require(auto_psf != benchmark.methods.end(),
          "Variable-seeing benchmark should estimate every frame PSF");
  require(with_psf.rgb_rmse < individual.rgb_rmse,
          "PSF-aware joint solve should beat individual demosaic RGB error");
  require(with_psf.star_fwhm_relative_error < without_psf.star_fwhm_relative_error,
          "PSF-aware forward model should improve FWHM recovery");
  require(with_psf.chroma_mae < without_psf.chroma_mae,
          "PSF-aware forward model should improve chroma recovery");
  require(auto_psf->metrics.star_fwhm_relative_error <
              without_psf.star_fwhm_relative_error,
          "Auto-estimated PSF should improve FWHM over the PSF-blind solve");
}

} // namespace

int main() {
  try {
    synthetic_scene_has_truth_and_cfa();
    benchmark_metrics_are_finite_and_accountable();
    star_chroma_guard_reduces_synthetic_star_false_color();
    morphological_poc_is_sensor_domain_and_finite();
    shared_profile_improves_non_gaussian_psf_shape();
    shared_epsf_recovers_anisotropic_psf();
    chromatic_epsf_reduces_lateral_color_error();
    multiframe_benchmark_is_reproducible_and_auditable();
    psf_aware_solver_recovers_variable_seeing_detail();
  } catch(const std::exception &error) {
    std::cerr << "synthetic_benchmark_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
