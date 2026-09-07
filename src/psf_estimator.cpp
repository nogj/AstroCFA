#include "astrocfa/psf_estimator.hpp"

#include "astrocfa/luminance_proxy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

double vector_median(std::vector<double> values) {
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

bool local_maximum(const astrocfa::LuminanceProxy &proxy, std::size_t x,
                   std::size_t y) {
  const float center = proxy.value(x, y);
  for(int dy = -1; dy <= 1; ++dy) {
    for(int dx = -1; dx <= 1; ++dx) {
      if(dx == 0 && dy == 0) {
        continue;
      }
      const std::size_t nx = static_cast<std::size_t>(static_cast<int>(x) + dx);
      const std::size_t ny = static_cast<std::size_t>(static_cast<int>(y) + dy);
      if(proxy.value(nx, ny) >= center) {
        return false;
      }
    }
  }
  return true;
}

int channel_index(astrocfa::CfaColor color) {
  if(color == astrocfa::CfaColor::red) {
    return 0;
  }
  if(color == astrocfa::CfaColor::blue) {
    return 2;
  }
  return 1;
}

struct LinearFitSums {
  double count = 0.0;
  double profile = 0.0;
  double value = 0.0;
  double profile_square = 0.0;
  double profile_value = 0.0;
  double value_square = 0.0;
};

double fit_error(const LinearFitSums &sums) {
  const double denominator =
      sums.count * sums.profile_square - sums.profile * sums.profile;
  if(sums.count < 4.0 || denominator <= 1.0e-12) {
    return std::numeric_limits<double>::infinity();
  }
  const double amplitude =
      (sums.count * sums.profile_value - sums.profile * sums.value) /
      denominator;
  if(amplitude <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  const double background =
      (sums.value - amplitude * sums.profile) / sums.count;
  return sums.value_square + amplitude * amplitude * sums.profile_square +
         sums.count * background * background +
         2.0 * amplitude * background * sums.profile -
         2.0 * amplitude * sums.profile_value -
         2.0 * background * sums.value;
}

double fit_cfa_gaussian_sigma(const astrocfa::CfaFrame &cfa, double initial_x,
                              double initial_y, std::size_t fit_radius) {
  const int radius = static_cast<int>(2U * fit_radius);
  const int center_x = static_cast<int>(std::lround(initial_x));
  const int center_y = static_cast<int>(std::lround(initial_y));
  double best_sigma = 0.0;
  double best_error = std::numeric_limits<double>::infinity();
  for(int offset_y = -4; offset_y <= 4; ++offset_y) {
    for(int offset_x = -4; offset_x <= 4; ++offset_x) {
      const double centroid_x = initial_x + 0.25 * offset_x;
      const double centroid_y = initial_y + 0.25 * offset_y;
      for(int sigma_step = 4; sigma_step <= 40; ++sigma_step) {
        const double sigma = 0.1 * sigma_step;
        LinearFitSums sums[3];
        for(int dy = -radius; dy <= radius; ++dy) {
          const int y = center_y + dy;
          if(y < 0 || y >= static_cast<int>(cfa.height())) {
            continue;
          }
          for(int dx = -radius; dx <= radius; ++dx) {
            const int x = center_x + dx;
            if(x < 0 || x >= static_cast<int>(cfa.width())) {
              continue;
            }
            const astrocfa::CfaSample sample = cfa.sample_info(
                static_cast<std::size_t>(x), static_cast<std::size_t>(y));
            if(!sample.valid || sample.clipped) {
              continue;
            }
            const double centered_x = static_cast<double>(x) - centroid_x;
            const double centered_y = static_cast<double>(y) - centroid_y;
            const double profile = std::exp(
                -0.5 * (centered_x * centered_x + centered_y * centered_y) /
                (sigma * sigma));
            LinearFitSums &fit =
                sums[channel_index(cfa.pattern().at(
                    static_cast<std::size_t>(x), static_cast<std::size_t>(y)))];
            fit.count += 1.0;
            fit.profile += profile;
            fit.value += sample.value;
            fit.profile_square += profile * profile;
            fit.profile_value += profile * sample.value;
            fit.value_square += sample.value * sample.value;
          }
        }
        double error = 0.0;
        double count = 0.0;
        for(const LinearFitSums &fit : sums) {
          const double channel_error = fit_error(fit);
          if(!std::isfinite(channel_error)) {
            error = std::numeric_limits<double>::infinity();
            break;
          }
          error += std::max(0.0, channel_error);
          count += fit.count;
        }
        if(count > 0.0 && error / count < best_error) {
          best_error = error / count;
          best_sigma = sigma;
        }
      }
    }
  }
  return best_sigma;
}

astrocfa::PsfStarMeasurement measure_star(
    const astrocfa::LuminanceProxy &proxy, std::size_t peak_x,
    std::size_t peak_y, double fallback_background, double background_sigma,
    const astrocfa::PsfEstimationOptions &options) {
  const int radius = static_cast<int>(options.fit_radius);
  const int annulus_outer = radius + 3;
  std::vector<double> annulus;
  for(int dy = -annulus_outer; dy <= annulus_outer; ++dy) {
    for(int dx = -annulus_outer; dx <= annulus_outer; ++dx) {
      const int radius_squared = dx * dx + dy * dy;
      if(radius_squared < (radius + 1) * (radius + 1) ||
         radius_squared > annulus_outer * annulus_outer) {
        continue;
      }
      const std::size_t x =
          static_cast<std::size_t>(static_cast<int>(peak_x) + dx);
      const std::size_t y =
          static_cast<std::size_t>(static_cast<int>(peak_y) + dy);
      annulus.push_back(proxy.value(x, y));
    }
  }
  const double background =
      annulus.empty() ? fallback_background : vector_median(std::move(annulus));
  double flux = 0.0;
  double weighted_x = 0.0;
  double weighted_y = 0.0;
  for(int dy = -radius; dy <= radius; ++dy) {
    for(int dx = -radius; dx <= radius; ++dx) {
      if(dx * dx + dy * dy > radius * radius) {
        continue;
      }
      const std::size_t x =
          static_cast<std::size_t>(static_cast<int>(peak_x) + dx);
      const std::size_t y =
          static_cast<std::size_t>(static_cast<int>(peak_y) + dy);
      const double signal = std::max(0.0, proxy.value(x, y) - background);
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
  double xx = 0.0;
  double yy = 0.0;
  double xy = 0.0;
  for(int dy = -radius; dy <= radius; ++dy) {
    for(int dx = -radius; dx <= radius; ++dx) {
      if(dx * dx + dy * dy > radius * radius) {
        continue;
      }
      const std::size_t x =
          static_cast<std::size_t>(static_cast<int>(peak_x) + dx);
      const std::size_t y =
          static_cast<std::size_t>(static_cast<int>(peak_y) + dy);
      const double signal = std::max(0.0, proxy.value(x, y) - background);
      const double centered_x = static_cast<double>(x) - centroid_x;
      const double centered_y = static_cast<double>(y) - centroid_y;
      xx += signal * centered_x * centered_x;
      yy += signal * centered_y * centered_y;
      xy += signal * centered_x * centered_y;
    }
  }
  xx /= flux;
  yy /= flux;
  xy /= flux;
  const double trace = xx + yy;
  const double discriminant = std::sqrt(
      std::max(0.0, (xx - yy) * (xx - yy) + 4.0 * xy * xy));
  const double major = 0.5 * (trace + discriminant);
  const double minor = 0.5 * (trace - discriminant);
  if(minor <= 1.0e-6 || major <= 0.0) {
    return {};
  }

  // Convert half-resolution proxy moments to sensor pixels and remove the
  // variance of the 2-pixel box integration used to build the proxy.
  const double sensor_variance = std::max(0.04, 2.0 * trace - 1.0 / 3.0);
  const double sigma = std::sqrt(sensor_variance);
  const double peak = proxy.value(peak_x, peak_y) - background;
  return astrocfa::PsfStarMeasurement{
      .x = 2.0 * (centroid_x + 0.5),
      .y = 2.0 * (centroid_y + 0.5),
      .sigma = sigma,
      .fwhm = 2.354820045 * sigma,
      .elongation = std::sqrt(major / minor),
      .snr = peak / std::max(background_sigma, 1.0e-6),
  };
}

} // namespace

namespace astrocfa {

PsfEstimate estimate_cfa_psf(const CfaFrame &cfa,
                             PsfEstimationOptions options) {
  if(!std::isfinite(options.detection_sigma) ||
     !std::isfinite(options.minimum_snr) ||
     !std::isfinite(options.maximum_elongation) ||
     !std::isfinite(options.maximum_peak) || options.detection_sigma <= 0.0 ||
     options.minimum_snr <= 0.0 || options.maximum_elongation < 1.0 ||
     options.maximum_peak <= 0.0 || options.maximum_peak > 1.25 ||
     options.fit_radius < 2 ||
     options.maximum_stars == 0) {
    throw std::invalid_argument("Invalid CFA PSF estimation options");
  }
  const LuminanceProxy proxy = build_cfa_luminance_proxy(cfa);
  PsfEstimate estimate;
  const std::size_t margin = options.fit_radius + 3U;
  if(proxy.width() <= 2U * margin || proxy.height() <= 2U * margin) {
    return estimate;
  }

  std::vector<float> samples;
  samples.reserve(proxy.width() * proxy.height());
  for(std::size_t y = 0; y < proxy.height(); ++y) {
    for(std::size_t x = 0; x < proxy.width(); ++x) {
      const float value = proxy.value(x, y);
      if(value > 0.0F) {
        samples.push_back(value);
      }
    }
  }
  if(samples.empty()) {
    return estimate;
  }
  const RobustBackground background = robust_background(samples);
  const double threshold =
      background.location + options.detection_sigma * background.sigma;

  struct Candidate {
    std::size_t x = 0;
    std::size_t y = 0;
    float peak = 0.0F;
  };
  std::vector<Candidate> candidates;
  for(std::size_t y = margin; y + margin < proxy.height(); ++y) {
    for(std::size_t x = margin; x + margin < proxy.width(); ++x) {
      const float peak = proxy.value(x, y);
      if(peak > threshold && peak < options.maximum_peak &&
         local_maximum(proxy, x, y)) {
        candidates.push_back(Candidate{.x = x, .y = y, .peak = peak});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &left, const Candidate &right) {
              return left.peak > right.peak;
            });
  estimate.candidates = candidates.size();

  for(const Candidate &candidate : candidates) {
    if(estimate.stars.size() >= options.maximum_stars) {
      break;
    }
    const PsfStarMeasurement star =
        measure_star(proxy, candidate.x, candidate.y, background.location,
                     background.sigma, options);
    PsfStarMeasurement fitted_star = star;
    fitted_star.sigma =
        fit_cfa_gaussian_sigma(cfa, star.x, star.y, options.fit_radius);
    fitted_star.fwhm = 2.354820045 * fitted_star.sigma;
    if(fitted_star.sigma <= 0.0 || fitted_star.snr < options.minimum_snr ||
       fitted_star.elongation > options.maximum_elongation) {
      continue;
    }
    bool isolated = true;
    for(const PsfStarMeasurement &accepted : estimate.stars) {
      if(std::hypot(fitted_star.x - accepted.x, fitted_star.y - accepted.y) <
         4.0 * static_cast<double>(options.fit_radius)) {
        isolated = false;
        break;
      }
    }
    if(isolated) {
      estimate.stars.push_back(fitted_star);
    }
  }
  if(estimate.stars.empty()) {
    return estimate;
  }

  std::vector<double> sigmas;
  std::vector<double> elongations;
  sigmas.reserve(estimate.stars.size());
  elongations.reserve(estimate.stars.size());
  for(const PsfStarMeasurement &star : estimate.stars) {
    sigmas.push_back(star.sigma);
    elongations.push_back(star.elongation);
  }
  const double center = vector_median(sigmas);
  std::vector<double> deviations;
  deviations.reserve(sigmas.size());
  for(double sigma : sigmas) {
    deviations.push_back(std::abs(sigma - center));
  }
  const double scatter = 1.4826 * vector_median(deviations);
  std::vector<double> inliers;
  for(double sigma : sigmas) {
    if(scatter <= 1.0e-9 || std::abs(sigma - center) <= 3.0 * scatter) {
      inliers.push_back(sigma);
    }
  }
  estimate.sigma = vector_median(inliers);
  estimate.fwhm = 2.354820045 * estimate.sigma;
  estimate.median_elongation = vector_median(elongations);
  estimate.scatter = scatter;
  estimate.used_stars = inliers.size();
  estimate.valid = estimate.used_stars > 0;
  return estimate;
}

std::vector<double> relative_psf_sigmas(
    const std::vector<PsfEstimate> &estimates, double strength) {
  if(strength < 0.0 || strength > 1.0 || !std::isfinite(strength)) {
    throw std::invalid_argument("Relative PSF strength must be in 0..1");
  }
  double minimum_variance = std::numeric_limits<double>::infinity();
  for(const PsfEstimate &estimate : estimates) {
    if(!estimate.valid || estimate.sigma <= 0.0) {
      throw std::invalid_argument(
          "Relative PSF requires a valid estimate per frame");
    }
    minimum_variance =
        std::min(minimum_variance, estimate.sigma * estimate.sigma);
  }
  std::vector<double> relative;
  relative.reserve(estimates.size());
  for(const PsfEstimate &estimate : estimates) {
    relative.push_back(
        strength * std::sqrt(std::max(
                       0.0, estimate.sigma * estimate.sigma - minimum_variance)));
  }
  return relative;
}

} // namespace astrocfa
