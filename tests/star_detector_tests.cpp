#include "astrocfa/psf_estimator.hpp"
#include "astrocfa/star_detector.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void detects_bright_cfa_safe_star_candidate() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(8, 8, rggb);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(
          x, y,
          astrocfa::CfaSample{.value = 0.01F, .valid = true, .clipped = false});
    }
  }

  for(std::size_t y = 2; y < 4; ++y) {
    for(std::size_t x = 2; x < 4; ++x) {
      frame.set_sample(x, y, astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = false});
    }
  }

  const astrocfa::StarDetectionStats stats = astrocfa::detect_star_candidates(frame, 2.0);
  require(stats.candidates == 1, "One bright synthetic star should be detected");
  require(stats.largest_area == 1, "The compact 2x2 star maps to one proxy sample");
  require(stats.brightest > 0.9F, "Brightest proxy sample should preserve stellar signal");
}

void ignores_clipped_blocks_when_estimating_background() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame frame(4, 4, rggb);
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      frame.set_sample(
          x, y,
          astrocfa::CfaSample{.value = 0.1F, .valid = true, .clipped = false});
    }
  }

  frame.set_sample(
      0, 0,
      astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(
      1, 0,
      astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(
      0, 1,
      astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});
  frame.set_sample(
      1, 1,
      astrocfa::CfaSample{.value = 1.0F, .valid = true, .clipped = true});

  const astrocfa::StarDetectionStats stats = astrocfa::detect_star_candidates(frame, 3.0);
  require(stats.background_mean < 0.11,
          "Clipped block should not bias background upward");
}

astrocfa::CfaFrame gaussian_star_field(double sigma) {
  astrocfa::CfaFrame frame(64, 48, astrocfa::BayerPattern{});
  const double stars[][2] = {{14.0, 12.0}, {34.0, 14.0},
                             {20.0, 34.0}, {48.0, 32.0}};
  for(std::size_t y = 0; y < frame.height(); ++y) {
    for(std::size_t x = 0; x < frame.width(); ++x) {
      double value = 0.01;
      for(const auto &star : stars) {
        const double dx = static_cast<double>(x) - star[0];
        const double dy = static_cast<double>(y) - star[1];
        value += 0.55 * std::exp(-0.5 * (dx * dx + dy * dy) /
                                 (sigma * sigma));
      }
      frame.set_sample(x, y, static_cast<float>(value));
    }
  }
  return frame;
}

void estimates_psf_from_cfa_safe_star_profiles() {
  const astrocfa::PsfEstimate estimate =
      astrocfa::estimate_cfa_psf(gaussian_star_field(2.0));
  require(estimate.valid, "Gaussian stars should produce a valid PSF estimate");
  require(estimate.used_stars >= 3, "PSF estimate should combine several stars");
  require(std::abs(estimate.sigma - 2.0) < 0.45,
          "Estimated PSF sigma should track the synthetic profile");
  require(estimate.median_elongation < 1.15,
          "Circular synthetic stars should remain nearly circular");
}

void converts_absolute_psfs_to_sharpest_relative_kernel() {
  const std::vector<double> relative = astrocfa::relative_psf_sigmas(
      {{.valid = true, .sigma = 2.0}, {.valid = true, .sigma = 2.5}}, 1.0);
  require(relative.size() == 2, "Relative PSF should preserve frame count");
  require(relative[0] == 0.0, "Sharpest frame should define the latent PSF");
  require(std::abs(relative[1] - 1.5) < 1.0e-12,
          "Relative Gaussian variance should subtract in quadrature");
}

} // namespace

int main() {
  try {
    detects_bright_cfa_safe_star_candidate();
    ignores_clipped_blocks_when_estimating_background();
    estimates_psf_from_cfa_safe_star_profiles();
    converts_absolute_psfs_to_sharpest_relative_kernel();
  } catch(const std::exception &error) {
    std::cerr << "star_detector_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
