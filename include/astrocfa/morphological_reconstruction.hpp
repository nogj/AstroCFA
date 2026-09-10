#pragma once

#include "astrocfa/demosaic.hpp"

#include <cstddef>
#include <limits>

namespace astrocfa {

struct MorphologicalReconstructionOptions {
  double diffuse_sigma = 4.0;
  double detection_sigma = 6.0;
  double source_mask_radius = 4.0;
  double minimum_source_sigma = 0.35;
  double maximum_source_sigma = 1.30;
  double source_sigma_step = 0.10;
  double position_search_radius = 1.0;
  double position_search_step = 0.125;
  std::size_t maximum_sources = 64;
  bool share_psf_across_sources = true;
  bool enable_epsf = true;
  bool enable_chromatic_epsf = true;
  double profile_bin_size = 0.125;
  double profile_radius = 4.0;
  std::size_t profile_iterations = 4;
  std::size_t maximum_profile_sources = 12;
  std::size_t epsf_oversampling = 2;
  std::size_t epsf_iterations = 8;
};

enum class MorphologicalModel {
  independent_gaussian,
  shared_gaussian,
  shared_profile,
  shared_epsf,
  shared_chromatic_epsf,
};

struct MorphologicalReconstructionStats {
  std::size_t detected_sources = 0;
  std::size_t validated_sources = 0;
  std::size_t reconstructed_sources = 0;
  std::size_t profile_sources = 0;
  MorphologicalModel selected_model = MorphologicalModel::independent_gaussian;
  double independent_validation = std::numeric_limits<double>::infinity();
  double shared_validation = std::numeric_limits<double>::infinity();
  double profile_validation = std::numeric_limits<double>::infinity();
  double epsf_validation = std::numeric_limits<double>::infinity();
  double chromatic_epsf_validation = std::numeric_limits<double>::infinity();
  double chromatic_complexity_penalty = 0.0;
  double profile_information_gain = 0.0;
  double red_epsf_shift_x = 0.0;
  double red_epsf_shift_y = 0.0;
  double red_epsf_scale = 1.0;
  double blue_epsf_shift_x = 0.0;
  double blue_epsf_shift_y = 0.0;
  double blue_epsf_scale = 1.0;
};

struct MorphologicalReconstructionResult {
  DemosaicResult reconstruction;
  MorphologicalReconstructionStats stats;
};

[[nodiscard]] const char *morphological_model_name(MorphologicalModel model);

[[nodiscard]] MorphologicalReconstructionResult
reconstruct_morphological_cfa_detailed(
    const CfaFrame &cfa, const NoiseModel &noise_model,
    MorphologicalReconstructionOptions options = {});

// Direct CFA inversion with separately modeled diffuse and point-source components.
[[nodiscard]] DemosaicResult reconstruct_morphological_cfa(
    const CfaFrame &cfa, const NoiseModel &noise_model,
    MorphologicalReconstructionOptions options = {});

} // namespace astrocfa
