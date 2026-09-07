#include "astrocfa/multiframe_benchmark.hpp"

#include "astrocfa/demosaic.hpp"
#include "astrocfa/synthetic_astro_scene.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>

namespace {

astrocfa::RgbPixel sample_rgb(const astrocfa::RgbImage &image, double x, double y,
                              bool &valid) {
  valid = x >= 0.0 && y >= 0.0 && x <= static_cast<double>(image.width() - 1U) &&
          y <= static_cast<double>(image.height() - 1U);
  if(!valid) {
    return {};
  }
  const std::size_t x0 = static_cast<std::size_t>(x);
  const std::size_t y0 = static_cast<std::size_t>(y);
  const std::size_t x1 = std::min(x0 + 1U, image.width() - 1U);
  const std::size_t y1 = std::min(y0 + 1U, image.height() - 1U);
  const double fx = x - static_cast<double>(x0);
  const double fy = y - static_cast<double>(y0);
  const double weights[4] = {
      (1.0 - fx) * (1.0 - fy), fx * (1.0 - fy),
      (1.0 - fx) * fy, fx * fy,
  };
  const astrocfa::RgbPixel pixels[4] = {
      image.pixel(x0, y0), image.pixel(x1, y0),
      image.pixel(x0, y1), image.pixel(x1, y1),
  };
  astrocfa::RgbPixel result;
  for(int i = 0; i < 4; ++i) {
    result.r += static_cast<float>(weights[i] * pixels[i].r);
    result.g += static_cast<float>(weights[i] * pixels[i].g);
    result.b += static_cast<float>(weights[i] * pixels[i].b);
  }
  return result;
}

astrocfa::RgbImage aligned_individual_average(
    const std::vector<astrocfa::SyntheticCfaObservation> &observations,
    const std::vector<astrocfa::SubpixelOffset> &offsets) {
  std::vector<astrocfa::RgbImage> reconstructions;
  reconstructions.reserve(observations.size());
  for(const auto &observation : observations) {
    reconstructions.push_back(astrocfa::reconstruct_inverse_refine(
                                  observation.cfa, astrocfa::NoiseModel{},
                                  astrocfa::InverseRefinementOptions{
                                      .frequency = {.tile_size = 16},
                                      .iterations = 3,
                                  })
                                  .image);
  }

  astrocfa::RgbImage average(reconstructions.front().width(),
                            reconstructions.front().height());
  for(std::size_t y = 0; y < average.height(); ++y) {
    for(std::size_t x = 0; x < average.width(); ++x) {
      astrocfa::RgbPixel sum;
      std::size_t count = 0;
      for(std::size_t i = 0; i < reconstructions.size(); ++i) {
        bool valid = false;
        const astrocfa::RgbPixel sample = sample_rgb(
            reconstructions[i], static_cast<double>(x) - offsets[i].dx,
            static_cast<double>(y) - offsets[i].dy, valid);
        if(valid) {
          sum.r += sample.r;
          sum.g += sample.g;
          sum.b += sample.b;
          count += 1;
        }
      }
      if(count > 0) {
        const float divisor = static_cast<float>(count);
        sum.r /= divisor;
        sum.g /= divisor;
        sum.b /= divisor;
      }
      average.set_pixel(x, y, sum);
    }
  }
  return average;
}

std::vector<astrocfa::SubpixelOffset> benchmark_offsets(std::size_t count) {
  const astrocfa::SubpixelOffset sequence[] = {
      {.dx = 0.0, .dy = 0.0},
      {.dx = 0.5, .dy = 0.0},
      {.dx = 0.0, .dy = 0.5},
      {.dx = 0.5, .dy = 0.5},
      {.dx = 0.25, .dy = 0.75},
      {.dx = 0.75, .dy = 0.25},
      {.dx = -0.25, .dy = 0.5},
      {.dx = 0.5, .dy = -0.25},
  };
  std::vector<astrocfa::SubpixelOffset> offsets;
  offsets.reserve(count);
  for(std::size_t i = 0; i < count; ++i) {
    offsets.push_back(sequence[i % (sizeof(sequence) / sizeof(sequence[0]))]);
  }
  return offsets;
}

} // namespace

namespace astrocfa {

MultiframeBenchmarkResult run_multiframe_benchmark(MultiframeBenchmarkOptions options) {
  if(options.width < 16 || options.height < 16 || options.frames == 0 ||
     options.frames > 32) {
    throw std::invalid_argument("Multiframe benchmark requires >=16x16 and 1..32 frames");
  }

  const SyntheticAstroScene scene = make_synthetic_astro_scene(
      SyntheticAstroSceneOptions{
          .width = options.width,
          .height = options.height,
          .seed = options.seed,
          .add_noise = false,
          .add_hot_pixels = false,
      });
  std::vector<SubpixelOffset> offsets = benchmark_offsets(options.frames);
  std::vector<SyntheticCfaObservation> observations;
  observations.reserve(options.frames);
  std::size_t injected_transients = 0;
  for(std::size_t i = 0; i < options.frames; ++i) {
    const double seeing = options.vary_seeing ? 0.35 + 0.12 * (i % 4U) : 0.0;
    observations.push_back(make_synthetic_cfa_observation(
        scene.truth, SyntheticObservationOptions{
                         .offset = offsets[i],
                         .psf_sigma = seeing,
                         .add_noise = options.add_noise,
                         .transient_samples = options.transients_per_frame,
                         .seed = options.seed + static_cast<std::uint32_t>(101U * i),
                     }));
    injected_transients += observations.back().transient_samples;
  }

  std::vector<JointCfaFrame> joint_inputs;
  joint_inputs.reserve(observations.size());
  for(std::size_t i = 0; i < observations.size(); ++i) {
    joint_inputs.push_back(JointCfaFrame{
        .cfa = &observations[i].cfa,
        .offset = offsets[i],
    });
  }

  MultiframeBenchmarkResult benchmark{
      .truth = scene.truth,
      .stars = scene.stars,
      .offsets = offsets,
      .injected_transients = injected_transients,
  };
  RgbImage individual = aligned_individual_average(observations, offsets);
  benchmark.methods.push_back(MultiframeBenchmarkMethod{
      .name = "individual-inverse-average",
      .image = std::move(individual),
  });

  const JointReconstructionOptions base_options{
      .iterations = 0,
      .luma_smoothness = options.luma_smoothness,
      .chroma_smoothness = options.chroma_smoothness,
      .noise = NoiseModel{.read_noise = 0.0025, .shot_noise_scale = 0.0018},
  };
  JointReconstructionResult initialization =
      reconstruct_joint_cfa(joint_inputs, base_options);
  benchmark.methods.push_back(MultiframeBenchmarkMethod{
      .name = "cfa-splat-initialization",
      .image = std::move(initialization.image),
      .confidence = std::make_unique<RgbImage>(std::move(initialization.confidence)),
      .solver_stats = initialization.stats,
      .has_solver_stats = true,
  });

  JointReconstructionOptions least_squares_options = base_options;
  least_squares_options.iterations = options.iterations;
  least_squares_options.huber_sigma = 1.0e9;
  JointReconstructionResult least_squares =
      reconstruct_joint_cfa(joint_inputs, least_squares_options);
  benchmark.methods.push_back(MultiframeBenchmarkMethod{
      .name = "joint-no-robust",
      .image = std::move(least_squares.image),
      .confidence = std::make_unique<RgbImage>(std::move(least_squares.confidence)),
      .solver_stats = least_squares.stats,
      .has_solver_stats = true,
  });

  JointReconstructionOptions robust_options = base_options;
  robust_options.iterations = options.iterations;
  JointReconstructionResult robust = reconstruct_joint_cfa(joint_inputs, robust_options);
  benchmark.methods.push_back(MultiframeBenchmarkMethod{
      .name = "joint-robust",
      .image = std::move(robust.image),
      .confidence = std::make_unique<RgbImage>(std::move(robust.confidence)),
      .solver_stats = robust.stats,
      .has_solver_stats = true,
  });

  for(auto &method : benchmark.methods) {
    method.metrics = measure_reconstruction(scene.truth, method.image,
                                            observations.front().cfa, scene.stars);
  }
  return benchmark;
}

} // namespace astrocfa
