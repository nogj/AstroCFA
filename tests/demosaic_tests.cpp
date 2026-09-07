#include "astrocfa/demosaic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void baseline_preserves_measured_cfa_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(4, 4, rggb);

  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = static_cast<float>(10U * y + x),
                               .valid = true,
                               .clipped = false,
                           });
    }
  }

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_baseline(cfa, astrocfa::NoiseModel{});
  require(result.residual.samples == 16, "Baseline should compare every valid sample");
  require(result.residual.maximum_absolute == 0.0F,
          "Baseline demosaic must preserve measured CFA samples exactly");
  require(result.noise_weighted_residual.reduced_chi_square == 0.0,
          "Exact remosaic should have zero weighted residual");
}

void baseline_skips_clipped_residual_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(2, 2, rggb);
  cfa.set_sample(0, 0, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  cfa.set_sample(1, 0, astrocfa::CfaSample{.value = 0.2F, .valid = true, .clipped = false});
  cfa.set_sample(0, 1, astrocfa::CfaSample{.value = 0.2F, .valid = true, .clipped = false});
  cfa.set_sample(1, 1, astrocfa::CfaSample{.value = 0.2F, .valid = true, .clipped = false});

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_baseline(cfa, astrocfa::NoiseModel{});
  require(result.residual.samples == 3, "Clipped samples should be excluded from residual");
}

void frequency_guided_preserves_measured_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const float value = ((x + y) % 2U) == 0U ? 1.0F : 0.0F;
      cfa.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true, .clipped = false});
    }
  }

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_frequency_guided(cfa, astrocfa::NoiseModel{},
                                             astrocfa::FrequencyGuidedDemosaicOptions{
                                                 .frequency = {.tile_size = 4},
                                                 .max_chroma_suppression = 1.0,
                                             });
  require(result.residual.maximum_absolute == 0.0F,
          "Frequency-guided demosaic must preserve measured CFA samples");
}

void malvar_preserves_measured_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = static_cast<float>((x + y) % 5U) / 5.0F,
                               .valid = true,
                               .clipped = false,
                           });
    }
  }

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_malvar_baseline(cfa, astrocfa::NoiseModel{});
  require(result.residual.maximum_absolute == 0.0F,
          "MHC-like baseline must preserve measured CFA samples");
}

void residual_interpolation_preserves_measured_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = static_cast<float>((3U * x + 5U * y) % 11U) / 11.0F,
                               .valid = true,
                               .clipped = false,
                           });
    }
  }

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_residual_interpolation(cfa, astrocfa::NoiseModel{});
  require(result.residual.maximum_absolute == 0.0F,
          "Residual interpolation must preserve measured CFA samples");
}

void residual_interpolation_tracks_color_difference() {
  astrocfa::BayerPattern rggb;
  astrocfa::RgbImage truth(10, 10);
  for(std::size_t y = 0; y < truth.height(); ++y) {
    for(std::size_t x = 0; x < truth.width(); ++x) {
      const float g = 0.2F + static_cast<float>(x + y) * 0.01F;
      truth.set_pixel(x, y, astrocfa::RgbPixel{.r = g + 0.1F, .g = g, .b = g - 0.05F});
    }
  }

  const astrocfa::CfaFrame cfa = astrocfa::remosaic(truth, rggb);
  const astrocfa::RgbImage reconstructed = astrocfa::demosaic_residual_interpolation(cfa);
  const astrocfa::RgbPixel pixel = reconstructed.pixel(4, 4);
  require(std::abs((pixel.r - pixel.g) - 0.1F) < 0.08F,
          "Residual interpolation should approximately preserve smooth R-G residuals");
}

void frequency_guided_suppresses_interpolated_chroma() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(8, 8, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const float value = ((x + y) % 2U) == 0U ? 1.0F : 0.0F;
      cfa.set_sample(x, y, astrocfa::CfaSample{.value = value, .valid = true, .clipped = false});
    }
  }

  const astrocfa::RgbImage baseline = astrocfa::demosaic_bilinear_baseline(cfa);
  const astrocfa::RgbImage guided = astrocfa::demosaic_frequency_guided(
      cfa, astrocfa::FrequencyGuidedDemosaicOptions{
               .frequency = {.tile_size = 4},
               .max_chroma_suppression = 1.0,
           });

  const astrocfa::RgbPixel before = baseline.pixel(0, 0);
  const astrocfa::RgbPixel after = guided.pixel(0, 0);
  const double before_chroma = std::abs(before.r - before.g) + std::abs(before.b - before.g);
  const double after_chroma = std::abs(after.r - after.g) + std::abs(after.b - after.g);
  require(after_chroma < before_chroma,
          "High-risk tiles should suppress unsupported interpolated chroma");

  const astrocfa::DemosaicQuality baseline_quality =
      astrocfa::analyze_demosaic_quality(baseline, cfa);
  const astrocfa::DemosaicQuality guided_quality =
      astrocfa::analyze_demosaic_quality(guided, cfa);
  require(guided_quality.mean_interpolated_chroma <
              baseline_quality.mean_interpolated_chroma,
          "Frequency guidance should reduce interpolated chroma magnitude");
}

void inverse_refine_preserves_measured_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(10, 10, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = static_cast<float>((7U * x + 3U * y) % 17U) / 17.0F,
                               .valid = true,
                               .clipped = false,
                           });
    }
  }

  const astrocfa::DemosaicResult result =
      astrocfa::reconstruct_inverse_refine(cfa, astrocfa::NoiseModel{},
                                           astrocfa::InverseRefinementOptions{
                                               .frequency = {.tile_size = 4},
                                               .iterations = 4,
                                           });
  require(result.residual.maximum_absolute == 0.0F,
          "Inverse refinement must preserve measured CFA samples");
}

void inverse_refine_reduces_chroma_roughness_on_alias_pattern() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame cfa(12, 12, rggb);
  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const float base = 0.30F + static_cast<float>(x + y) * 0.005F;
      const float carrier = ((x + y) % 2U) == 0U ? 0.25F : -0.20F;
      cfa.set_sample(x, y, astrocfa::CfaSample{
                               .value = std::clamp(base + carrier, 0.0F, 1.0F),
                               .valid = true,
                               .clipped = false,
                           });
    }
  }

  const astrocfa::RgbImage guided = astrocfa::demosaic_frequency_guided(
      cfa, astrocfa::FrequencyGuidedDemosaicOptions{
               .frequency = {.tile_size = 4},
           });
  const astrocfa::RgbImage refined =
      astrocfa::demosaic_inverse_refine(cfa, astrocfa::InverseRefinementOptions{
                                                 .frequency = {.tile_size = 4},
                                                 .iterations = 6,
                                             });
  const astrocfa::DemosaicQuality guided_quality =
      astrocfa::analyze_demosaic_quality(guided, cfa);
  const astrocfa::DemosaicQuality refined_quality =
      astrocfa::analyze_demosaic_quality(refined, cfa);

  require(refined_quality.mean_chroma_roughness <=
              guided_quality.mean_chroma_roughness,
          "Inverse refinement should not increase chroma roughness on alias-heavy data");
}

} // namespace

int main() {
  try {
    baseline_preserves_measured_cfa_samples();
    baseline_skips_clipped_residual_samples();
    malvar_preserves_measured_samples();
    residual_interpolation_preserves_measured_samples();
    residual_interpolation_tracks_color_difference();
    frequency_guided_preserves_measured_samples();
    frequency_guided_suppresses_interpolated_chroma();
    inverse_refine_preserves_measured_samples();
    inverse_refine_reduces_chroma_roughness_on_alias_pattern();
  } catch(const std::exception &error) {
    std::cerr << "demosaic_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
