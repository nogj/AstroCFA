#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>
#include <vector>

namespace astrocfa {

struct PsfEstimationOptions {
  double detection_sigma = 6.0;
  double minimum_snr = 12.0;
  double maximum_elongation = 1.8;
  double maximum_peak = 0.92;
  std::size_t fit_radius = 3;
  std::size_t maximum_stars = 64;
};

struct PsfStarMeasurement {
  double x = 0.0;
  double y = 0.0;
  double sigma = 0.0;
  double fwhm = 0.0;
  double elongation = 0.0;
  double snr = 0.0;
};

struct PsfEstimate {
  bool valid = false;
  double sigma = 0.0;
  double fwhm = 0.0;
  double median_elongation = 0.0;
  double scatter = 0.0;
  std::size_t candidates = 0;
  std::size_t used_stars = 0;
  std::vector<PsfStarMeasurement> stars;
};

[[nodiscard]] PsfEstimate estimate_cfa_psf(
    const CfaFrame &cfa, PsfEstimationOptions options = {});

[[nodiscard]] std::vector<double> relative_psf_sigmas(
    const std::vector<PsfEstimate> &estimates, double strength = 0.75);

} // namespace astrocfa
