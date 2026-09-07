#include "astrocfa/cfa.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void require_near(double actual, double expected, double tolerance, const char *message) {
  if(std::abs(actual - expected) > tolerance) {
    std::cerr << message << ": actual=" << actual << " expected=" << expected << "\n";
    throw std::runtime_error(message);
  }
}

void remosaic_uses_bayer_phase() {
  astrocfa::BayerPattern rggb;
  astrocfa::RgbImage rgb(2, 2);
  rgb.set_pixel(0, 0, {10.0F, 20.0F, 30.0F});
  rgb.set_pixel(1, 0, {11.0F, 21.0F, 31.0F});
  rgb.set_pixel(0, 1, {12.0F, 22.0F, 32.0F});
  rgb.set_pixel(1, 1, {13.0F, 23.0F, 33.0F});

  const astrocfa::CfaFrame cfa = astrocfa::remosaic(rgb, rggb);
  require_near(cfa.sample(0, 0), 10.0, 0.0, "R phase should use red channel");
  require_near(cfa.sample(1, 0), 21.0, 0.0, "G1 phase should use green channel");
  require_near(cfa.sample(0, 1), 22.0, 0.0, "G2 phase should use green channel");
  require_near(cfa.sample(1, 1), 33.0, 0.0, "B phase should use blue channel");
}

void exact_reconstruction_has_zero_residual() {
  astrocfa::BayerPattern rggb;
  astrocfa::RgbImage rgb(3, 3);
  for(std::size_t y = 0; y < rgb.height(); ++y) {
    for(std::size_t x = 0; x < rgb.width(); ++x) {
      const float base = static_cast<float>(10U * y + x);
      rgb.set_pixel(x, y, {base + 1.0F, base + 2.0F, base + 3.0F});
    }
  }

  const astrocfa::CfaFrame measured = astrocfa::remosaic(rgb, rggb);
  const astrocfa::RemosaicResidual residual =
      astrocfa::compute_remosaic_residual(measured, rgb);

  require(residual.samples == 9, "Residual should count every CFA sample");
  require_near(residual.mean_absolute, 0.0, 0.0, "Exact remosaic mean residual");
  require_near(residual.root_mean_square, 0.0, 0.0, "Exact remosaic RMS residual");
  require_near(residual.maximum_absolute, 0.0, 0.0, "Exact remosaic max residual");
}

void residual_detects_unsupported_rgb_change() {
  astrocfa::BayerPattern rggb;
  astrocfa::RgbImage truth(2, 2);
  truth.set_pixel(0, 0, {1.0F, 0.0F, 0.0F});
  truth.set_pixel(1, 0, {0.0F, 2.0F, 0.0F});
  truth.set_pixel(0, 1, {0.0F, 3.0F, 0.0F});
  truth.set_pixel(1, 1, {0.0F, 0.0F, 4.0F});
  const astrocfa::CfaFrame measured = astrocfa::remosaic(truth, rggb);

  astrocfa::RgbImage altered = truth;
  altered.set_pixel(1, 1, {0.0F, 0.0F, 10.0F});

  const astrocfa::RemosaicResidual residual =
      astrocfa::compute_remosaic_residual(measured, altered);

  require_near(residual.maximum_absolute, 6.0, 0.0, "Residual max should catch B change");
  require_near(residual.mean_absolute, 1.5, 0.0, "Residual mean should average CFA error");
  require_near(residual.root_mean_square, 3.0, 0.0, "Residual RMS should measure energy");
}

void residual_ignores_invalid_and_clipped_samples() {
  astrocfa::BayerPattern rggb;
  astrocfa::RgbImage truth(2, 2);
  truth.set_pixel(0, 0, {1.0F, 0.0F, 0.0F});
  truth.set_pixel(1, 0, {0.0F, 2.0F, 0.0F});
  truth.set_pixel(0, 1, {0.0F, 3.0F, 0.0F});
  truth.set_pixel(1, 1, {0.0F, 0.0F, 4.0F});

  astrocfa::CfaFrame measured = astrocfa::remosaic(truth, rggb);
  measured.set_valid(0, 0, false);
  measured.set_clipped(1, 1, true);

  astrocfa::RgbImage altered = truth;
  altered.set_pixel(0, 0, {100.0F, 0.0F, 0.0F});
  altered.set_pixel(1, 1, {0.0F, 0.0F, 100.0F});

  const astrocfa::RemosaicResidual residual =
      astrocfa::compute_remosaic_residual(measured, altered);

  require(residual.samples == 2, "Residual should skip invalid and clipped CFA samples");
  require_near(residual.mean_absolute, 0.0, 0.0, "Skipped samples should not add error");
}

void noise_weighted_residual_uses_sensor_noise_scale() {
  astrocfa::BayerPattern rggb;
  astrocfa::CfaFrame measured(1, 1, rggb);
  measured.set_sample(0, 0, astrocfa::CfaSample{.value = 0.25F, .valid = true, .clipped = false});

  astrocfa::RgbImage reconstructed(1, 1);
  reconstructed.set_pixel(0, 0, {0.26F, 0.0F, 0.0F});

  const astrocfa::NoiseModel low_noise{.read_noise = 0.001, .shot_noise_scale = 0.0001};
  const astrocfa::NoiseModel high_noise{.read_noise = 0.01, .shot_noise_scale = 0.01};

  const astrocfa::NoiseWeightedResidual low =
      astrocfa::compute_noise_weighted_remosaic_residual(measured, reconstructed, low_noise);
  const astrocfa::NoiseWeightedResidual high =
      astrocfa::compute_noise_weighted_remosaic_residual(measured, reconstructed, high_noise);

  require(low.samples == 1, "Weighted residual should count measured samples");
  require(low.reduced_chi_square > high.reduced_chi_square,
          "Same reconstruction error should matter less under a noisier model");
}

void quantization_preserves_metadata_and_sensor_codes() {
  astrocfa::CfaFrame input(2, 1, astrocfa::BayerPattern{});
  input.set_sample(0, 0,
                   astrocfa::CfaSample{.value = 0.3F, .valid = false, .clipped = false});
  input.set_sample(1, 0,
                   astrocfa::CfaSample{.value = 1.2F, .valid = true, .clipped = true});

  const astrocfa::CfaFrame quantized = astrocfa::quantize_cfa(input, 255U);
  require_near(quantized.sample(0, 0), 77.0 / 255.0, 1e-7,
               "CFA quantization should round to a sensor code");
  require_near(quantized.sample(1, 0), 1.0, 0.0,
               "CFA quantization should clamp to the maximum code");
  require(!quantized.sample_info(0, 0).valid,
          "CFA quantization should preserve validity");
  require(quantized.sample_info(1, 0).clipped,
          "CFA quantization should preserve clipping");
}

} // namespace

int main() {
  try {
    remosaic_uses_bayer_phase();
    exact_reconstruction_has_zero_residual();
    residual_detects_unsupported_rgb_change();
    residual_ignores_invalid_and_clipped_samples();
    noise_weighted_residual_uses_sensor_noise_scale();
    quantization_preserves_metadata_and_sensor_codes();
  } catch(const std::exception &error) {
    std::cerr << "cfa_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
