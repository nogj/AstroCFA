#include "astrocfa/reconstruction_metrics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

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

struct StarMeasurement {
  double flux = 0.0;
  double fwhm = 0.0;
  double elongation = 1.0;
  bool valid = false;
};

double median(std::vector<double> values) {
  if(values.empty()) {
    return 0.0;
  }
  const std::size_t middle = values.size() / 2U;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if(values.size() % 2U != 0U) {
    return upper;
  }
  std::nth_element(values.begin(), values.begin() + middle - 1U, values.end());
  return 0.5 * (upper + values[middle - 1U]);
}

StarMeasurement measure_star(const astrocfa::RgbImage &image,
                             const astrocfa::SyntheticStar &star) {
  const double aperture_radius = std::max(3.0, 3.5 * star.sigma);
  const double annulus_inner = aperture_radius + 1.0;
  const double annulus_outer = annulus_inner + 2.0;
  const int x_begin = std::max(0, static_cast<int>(std::floor(star.x - annulus_outer)));
  const int y_begin = std::max(0, static_cast<int>(std::floor(star.y - annulus_outer)));
  const int x_end = std::min(static_cast<int>(image.width()) - 1,
                             static_cast<int>(std::ceil(star.x + annulus_outer)));
  const int y_end = std::min(static_cast<int>(image.height()) - 1,
                             static_cast<int>(std::ceil(star.y + annulus_outer)));

  std::vector<double> annulus;
  for(int y = y_begin; y <= y_end; ++y) {
    for(int x = x_begin; x <= x_end; ++x) {
      const double radius = std::hypot(static_cast<double>(x) - star.x,
                                       static_cast<double>(y) - star.y);
      if(radius >= annulus_inner && radius <= annulus_outer) {
        annulus.push_back(luma(image.pixel(static_cast<std::size_t>(x),
                                            static_cast<std::size_t>(y))));
      }
    }
  }
  if(annulus.empty()) {
    return {};
  }
  const double background = median(std::move(annulus));

  double flux = 0.0;
  double weighted_x = 0.0;
  double weighted_y = 0.0;
  for(int y = y_begin; y <= y_end; ++y) {
    for(int x = x_begin; x <= x_end; ++x) {
      const double dx = static_cast<double>(x) - star.x;
      const double dy = static_cast<double>(y) - star.y;
      if(dx * dx + dy * dy > aperture_radius * aperture_radius) {
        continue;
      }
      const double signal = std::max(
          0.0, luma(image.pixel(static_cast<std::size_t>(x),
                                static_cast<std::size_t>(y))) - background);
      flux += signal;
      weighted_x += signal * static_cast<double>(x);
      weighted_y += signal * static_cast<double>(y);
    }
  }
  if(flux <= std::numeric_limits<double>::epsilon()) {
    return {};
  }

  const double centroid_x = weighted_x / flux;
  const double centroid_y = weighted_y / flux;
  double moment_xx = 0.0;
  double moment_yy = 0.0;
  double moment_xy = 0.0;
  for(int y = y_begin; y <= y_end; ++y) {
    for(int x = x_begin; x <= x_end; ++x) {
      const double star_dx = static_cast<double>(x) - star.x;
      const double star_dy = static_cast<double>(y) - star.y;
      if(star_dx * star_dx + star_dy * star_dy >
         aperture_radius * aperture_radius) {
        continue;
      }
      const double signal = std::max(
          0.0, luma(image.pixel(static_cast<std::size_t>(x),
                                static_cast<std::size_t>(y))) - background);
      const double dx = static_cast<double>(x) - centroid_x;
      const double dy = static_cast<double>(y) - centroid_y;
      moment_xx += signal * dx * dx;
      moment_yy += signal * dy * dy;
      moment_xy += signal * dx * dy;
    }
  }
  moment_xx /= flux;
  moment_yy /= flux;
  moment_xy /= flux;
  const double trace = moment_xx + moment_yy;
  const double discriminant = std::sqrt(
      std::max(0.0, (moment_xx - moment_yy) * (moment_xx - moment_yy) +
                        4.0 * moment_xy * moment_xy));
  const double major_variance = std::max(0.0, 0.5 * (trace + discriminant));
  const double minor_variance = std::max(0.0, 0.5 * (trace - discriminant));
  if(major_variance <= std::numeric_limits<double>::epsilon()) {
    return {};
  }
  const double circular_sigma = std::sqrt(0.5 * trace);
  const double elongation = minor_variance > 1.0e-12
                                ? std::sqrt(major_variance / minor_variance)
                                : std::numeric_limits<double>::infinity();
  return StarMeasurement{
      .flux = flux,
      .fwhm = 2.354820045 * circular_sigma,
      .elongation = elongation,
      .valid = std::isfinite(elongation),
  };
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
  for(const SyntheticStar &star : stars) {
    const StarMeasurement expected = measure_star(truth, star);
    const StarMeasurement actual = measure_star(reconstructed, star);
    if(!expected.valid || !actual.valid || expected.flux <= 0.0 ||
       expected.fwhm <= 0.0) {
      continue;
    }
    metrics.star_flux_relative_error +=
        std::abs(actual.flux - expected.flux) / expected.flux;
    metrics.star_flux_relative_bias +=
        (actual.flux - expected.flux) / expected.flux;
    metrics.star_fwhm_relative_error +=
        std::abs(actual.fwhm - expected.fwhm) / expected.fwhm;
    metrics.star_elongation_error +=
        std::abs(actual.elongation - expected.elongation);
    metrics.measured_stars += 1;
  }
  if(metrics.measured_stars > 0) {
    const double count = static_cast<double>(metrics.measured_stars);
    metrics.star_flux_relative_error /= count;
    metrics.star_flux_relative_bias /= count;
    metrics.star_fwhm_relative_error /= count;
    metrics.star_elongation_error /= count;
  }
  metrics.cfa_residual_mae =
      compute_remosaic_residual(measured, reconstructed).mean_absolute;

  return metrics;
}

} // namespace astrocfa
