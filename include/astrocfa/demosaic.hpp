#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/frequency_cfa.hpp"

namespace astrocfa {

enum class DemosaicMethod {
  bilinear_baseline,
  malvar_baseline,
  residual_interpolation,
  frequency_guided,
  inverse_refine,
};

struct FrequencyGuidedDemosaicOptions {
  FrequencyCfaOptions frequency;
  double max_chroma_suppression = 0.85;
};

struct InverseRefinementOptions {
  FrequencyCfaOptions frequency;
  int iterations = 6;
  double chroma_smoothness = 0.55;
  double edge_sensitivity = 28.0;
  double alias_suppression = 0.45;
};

struct DemosaicResult {
  RgbImage image;
  RemosaicResidual residual;
  NoiseWeightedResidual noise_weighted_residual;
};

struct DemosaicQuality {
  double mean_chroma_roughness = 0.0;
  double max_chroma_roughness = 0.0;
  double mean_interpolated_chroma = 0.0;
  double max_interpolated_chroma = 0.0;
  std::size_t samples = 0;
  std::size_t interpolated_chroma_samples = 0;
};

[[nodiscard]] RgbImage demosaic_bilinear_baseline(const CfaFrame &cfa);
[[nodiscard]] RgbImage demosaic_malvar_baseline(const CfaFrame &cfa);
[[nodiscard]] RgbImage demosaic_residual_interpolation(const CfaFrame &cfa);
[[nodiscard]] RgbImage
demosaic_frequency_guided(const CfaFrame &cfa,
                          FrequencyGuidedDemosaicOptions options = {});
[[nodiscard]] RgbImage demosaic_inverse_refine(const CfaFrame &cfa,
                                               InverseRefinementOptions options = {});
[[nodiscard]] DemosaicResult reconstruct_baseline(const CfaFrame &cfa,
                                                  const NoiseModel &noise_model);
[[nodiscard]] DemosaicResult reconstruct_malvar_baseline(const CfaFrame &cfa,
                                                         const NoiseModel &noise_model);
[[nodiscard]] DemosaicResult reconstruct_residual_interpolation(const CfaFrame &cfa,
                                                                const NoiseModel &noise_model);
[[nodiscard]] DemosaicResult
reconstruct_frequency_guided(const CfaFrame &cfa, const NoiseModel &noise_model,
                             FrequencyGuidedDemosaicOptions options = {});
[[nodiscard]] DemosaicResult
reconstruct_inverse_refine(const CfaFrame &cfa, const NoiseModel &noise_model,
                           InverseRefinementOptions options = {});
[[nodiscard]] DemosaicQuality analyze_demosaic_quality(const RgbImage &image);
[[nodiscard]] DemosaicQuality analyze_demosaic_quality(const RgbImage &image,
                                                       const CfaFrame &cfa);

} // namespace astrocfa
