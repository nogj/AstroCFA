#include "astrocfa/joint_reconstruction.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

float measured_value(astrocfa::CfaColor color) {
  if(color == astrocfa::CfaColor::red) {
    return 0.20F;
  }
  if(color == astrocfa::CfaColor::blue) {
    return 0.10F;
  }
  return 0.35F;
}

astrocfa::CfaFrame constant_cfa(std::size_t width = 10, std::size_t height = 10) {
  astrocfa::CfaFrame frame(width, height, astrocfa::BayerPattern{});
  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      frame.set_sample(x, y, measured_value(frame.pattern().at(x, y)));
    }
  }
  return frame;
}

void robust_solver_rejects_transient_sample() {
  astrocfa::CfaFrame first = constant_cfa();
  astrocfa::CfaFrame second = constant_cfa();
  astrocfa::CfaFrame transient = constant_cfa();
  transient.set_sample(4, 4, 0.90F);

  const std::vector<astrocfa::JointCfaFrame> frames = {
      {.cfa = &first},
      {.cfa = &second},
      {.cfa = &transient},
  };
  const astrocfa::JointReconstructionResult result = astrocfa::reconstruct_joint_cfa(
      frames, astrocfa::JointReconstructionOptions{
                  .iterations = 8,
                  .chroma_smoothness = 0.0,
                  .stop_on_discrepancy = false,
              });

  require(std::abs(result.image.pixel(4, 4).r - 0.20F) < 0.03F,
          "Robust joint solve should follow repeated measurements, not a transient");
  require(result.stats.robust_outliers >= 1,
          "Transient measurement should remain auditable as an outlier");
  require(result.stats.frames == 3, "Joint solve should report all frames");
}

void single_frame_preserves_every_measured_cfa_sample() {
  astrocfa::CfaFrame frame = constant_cfa();
  const astrocfa::JointReconstructionResult result = astrocfa::reconstruct_joint_cfa(
      {{.cfa = &frame}}, astrocfa::JointReconstructionOptions{.iterations = 3});
  require(result.stats.final_rmse < 1.0e-7,
          "Regularization must not move directly measured CFA channels");
}

void dithers_increase_high_resolution_channel_coverage() {
  astrocfa::CfaFrame frame = constant_cfa(8, 8);
  const astrocfa::JointReconstructionOptions options{
      .scale = 2,
      .iterations = 1,
      .chroma_smoothness = 0.0,
  };
  const astrocfa::JointReconstructionResult single = astrocfa::reconstruct_joint_cfa(
      {{{.cfa = &frame}}}, options);
  const astrocfa::JointReconstructionResult dithered = astrocfa::reconstruct_joint_cfa(
      {
          {.cfa = &frame, .offset = {.dx = 0.0, .dy = 0.0}},
          {.cfa = &frame, .offset = {.dx = 0.5, .dy = 0.0}},
          {.cfa = &frame, .offset = {.dx = 0.0, .dy = 0.5}},
          {.cfa = &frame, .offset = {.dx = 0.5, .dy = 0.5}},
      },
      options);

  require(dithered.stats.channel_coverage[0] > single.stats.channel_coverage[0],
          "Dithering should add direct red support on the high-resolution grid");
  require(dithered.stats.channel_coverage[2] > single.stats.channel_coverage[2],
          "Dithering should add direct blue support on the high-resolution grid");
}

void invalid_psf_sigma_is_rejected() {
  astrocfa::CfaFrame frame = constant_cfa();
  bool rejected = false;
  try {
    (void)astrocfa::reconstruct_joint_cfa({{.cfa = &frame, .psf_sigma = 8.1}});
  } catch(const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "Unbounded PSF support should be rejected");
}

void psf_usage_is_reported() {
  astrocfa::CfaFrame frame = constant_cfa();
  const astrocfa::JointReconstructionResult result =
      astrocfa::reconstruct_joint_cfa(
          {{.cfa = &frame, .psf_sigma = 0.6},
           {.cfa = &frame, .psf_sigma = 1.1}},
          astrocfa::JointReconstructionOptions{.iterations = 1});
  require(result.stats.psf_frames == 2, "Every PSF-aware frame should be reported");
  require(std::abs(result.stats.minimum_psf_sigma - 0.6) < 1.0e-12,
          "Minimum PSF sigma should be reported");
  require(std::abs(result.stats.maximum_psf_sigma - 1.1) < 1.0e-12,
          "Maximum PSF sigma should be reported");
}

void discrepancy_stops_at_noise_consistent_solution() {
  astrocfa::CfaFrame first = constant_cfa();
  astrocfa::CfaFrame second = constant_cfa();
  const astrocfa::JointReconstructionResult result =
      astrocfa::reconstruct_joint_cfa(
          {{.cfa = &first}, {.cfa = &second}},
          astrocfa::JointReconstructionOptions{
              .iterations = 10,
              .luma_smoothness = 0.0,
              .chroma_smoothness = 0.0,
              .stop_on_discrepancy = true,
              .minimum_iterations = 1,
          });
  require(result.stats.stopped_by_discrepancy,
          "Noise-consistent solution should trigger discrepancy stopping");
  require(result.stats.iterations == 1,
          "Discrepancy stopping should honor the minimum iteration count");
  require(result.stats.final_reduced_chi_square < 1.0e-10,
          "Exact repeated data should have negligible reduced chi-square");
}

} // namespace

int main() {
  try {
    robust_solver_rejects_transient_sample();
    single_frame_preserves_every_measured_cfa_sample();
    dithers_increase_high_resolution_channel_coverage();
    invalid_psf_sigma_is_rejected();
    psf_usage_is_reported();
    discrepancy_stops_at_noise_consistent_solution();
  } catch(const std::exception &error) {
    std::cerr << "joint_reconstruction_tests failed: " << error.what() << "\n";
    return 1;
  }
  return 0;
}
