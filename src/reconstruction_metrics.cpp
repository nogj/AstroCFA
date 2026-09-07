#include "astrocfa/reconstruction_metrics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

double luma(astrocfa::RgbPixel pixel) {
  return 0.2126 * pixel.r + 0.7152 * pixel.g + 0.0722 * pixel.b;
}

double false_color(astrocfa::RgbPixel truth, astrocfa::RgbPixel reconstructed) {
  const double truth_rg = static_cast<double>(truth.r) - truth.g;
  const double truth_bg = static_cast<double>(truth.b) - truth.g;
  const double recon_rg = static_cast<double>(reconstructed.r) - reconstructed.g;
  const double recon_bg = static_cast<double>(reconstructed.b) - reconstructed.g;
  return 0.5 * (std::abs(recon_rg - truth_rg) + std::abs(recon_bg - truth_bg));
}

bool near_star(std::size_t x, std::size_t y, const std::vector<astrocfa::SyntheticStar> &stars) {
  for(const astrocfa::SyntheticStar &star : stars) {
    const double radius = std::max(2.0, star.sigma * 3.0);
    const double dx = static_cast<double>(x) - star.x;
    const double dy = static_cast<double>(y) - star.y;
    if(dx * dx + dy * dy <= radius * radius) {
      return true;
    }
  }
  return false;
}

} // namespace

namespace astrocfa {

ReconstructionMetrics measure_reconstruction(const RgbImage &truth,
                                             const RgbImage &reconstructed,
                                             const CfaFrame &measured,
                                             const std::vector<SyntheticStar> &stars) {
  if(truth.width() != reconstructed.width() || truth.height() != reconstructed.height() ||
     truth.width() != measured.width() || truth.height() != measured.height()) {
    throw std::invalid_argument("Metric input dimensions differ");
  }

  ReconstructionMetrics metrics;
  double absolute_sum = 0.0;
  double square_sum = 0.0;
  double chroma_sum = 0.0;
  double star_luma_square_sum = 0.0;

  for(std::size_t y = 0; y < truth.height(); ++y) {
    for(std::size_t x = 0; x < truth.width(); ++x) {
      const RgbPixel expected = truth.pixel(x, y);
      const RgbPixel actual = reconstructed.pixel(x, y);
      const double errors[3] = {
          static_cast<double>(actual.r) - expected.r,
          static_cast<double>(actual.g) - expected.g,
          static_cast<double>(actual.b) - expected.b,
      };
      for(double error : errors) {
        const double absolute = std::abs(error);
        absolute_sum += absolute;
        square_sum += error * error;
        metrics.max_abs_error = std::max(metrics.max_abs_error, absolute);
        metrics.samples += 1;
      }
      chroma_sum += false_color(expected, actual);

      if(near_star(x, y, stars)) {
        metrics.star_false_color += false_color(expected, actual);
        const double luma_error = luma(actual) - luma(expected);
        star_luma_square_sum += luma_error * luma_error;
        metrics.star_samples += 1;
      }
    }
  }

  if(metrics.samples > 0) {
    metrics.rgb_mae = absolute_sum / static_cast<double>(metrics.samples);
    metrics.rgb_rmse = std::sqrt(square_sum / static_cast<double>(metrics.samples));
  }
  const double pixels = static_cast<double>(truth.width() * truth.height());
  if(pixels > 0.0) {
    metrics.chroma_mae = chroma_sum / pixels;
  }
  if(metrics.star_samples > 0) {
    metrics.star_false_color /= static_cast<double>(metrics.star_samples);
    metrics.star_luma_rmse =
        std::sqrt(star_luma_square_sum / static_cast<double>(metrics.star_samples));
  }
  metrics.cfa_residual_mae =
      compute_remosaic_residual(measured, reconstructed).mean_absolute;

  return metrics;
}

} // namespace astrocfa
