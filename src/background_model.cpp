#include "astrocfa/background_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr std::size_t maximum_terms = 6;

struct TileSample {
  double x = 0.0;
  double y = 0.0;
  std::array<double, 3> value = {0.0, 0.0, 0.0};
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

std::array<double, maximum_terms> polynomial_basis(double x, double y) {
  return {1.0, x, y, x * x, x * y, y * y};
}

std::size_t term_count(std::size_t degree) {
  return degree == 0 ? 1U : degree == 1 ? 3U : 6U;
}

std::array<double, maximum_terms> solve_channel(
    const std::vector<TileSample> &samples, const std::vector<double> &weights,
    std::size_t channel, std::size_t terms) {
  double system[maximum_terms][maximum_terms + 1]{};
  for(std::size_t sample_index = 0; sample_index < samples.size();
      ++sample_index) {
    const auto basis =
        polynomial_basis(samples[sample_index].x, samples[sample_index].y);
    const double weight = weights[sample_index];
    for(std::size_t row = 0; row < terms; ++row) {
      for(std::size_t column = 0; column < terms; ++column) {
        system[row][column] += weight * basis[row] * basis[column];
      }
      system[row][terms] +=
          weight * basis[row] * samples[sample_index].value[channel];
    }
  }
  for(std::size_t diagonal = 0; diagonal < terms; ++diagonal) {
    system[diagonal][diagonal] += 1.0e-12;
  }

  for(std::size_t pivot = 0; pivot < terms; ++pivot) {
    std::size_t best = pivot;
    for(std::size_t row = pivot + 1U; row < terms; ++row) {
      if(std::abs(system[row][pivot]) > std::abs(system[best][pivot])) {
        best = row;
      }
    }
    if(std::abs(system[best][pivot]) < 1.0e-14) {
      throw std::runtime_error("Background polynomial fit is singular");
    }
    if(best != pivot) {
      for(std::size_t column = pivot; column <= terms; ++column) {
        std::swap(system[pivot][column], system[best][column]);
      }
    }
    const double divisor = system[pivot][pivot];
    for(std::size_t column = pivot; column <= terms; ++column) {
      system[pivot][column] /= divisor;
    }
    for(std::size_t row = 0; row < terms; ++row) {
      if(row == pivot) {
        continue;
      }
      const double factor = system[row][pivot];
      for(std::size_t column = pivot; column <= terms; ++column) {
        system[row][column] -= factor * system[pivot][column];
      }
    }
  }

  std::array<double, maximum_terms> coefficients{};
  for(std::size_t row = 0; row < terms; ++row) {
    coefficients[row] = system[row][terms];
  }
  return coefficients;
}

double evaluate(const std::array<double, maximum_terms> &coefficients,
                std::size_t terms, double x, double y) {
  const auto basis = polynomial_basis(x, y);
  double value = 0.0;
  for(std::size_t term = 0; term < terms; ++term) {
    value += coefficients[term] * basis[term];
  }
  return value;
}

std::vector<TileSample> collect_tile_samples(
    const astrocfa::RgbImage &image,
    const astrocfa::BackgroundModelOptions &options) {
  struct Pixel {
    double luminance = 0.0;
    std::array<double, 3> value = {0.0, 0.0, 0.0};
  };
  std::vector<TileSample> samples;
  for(std::size_t top = 0; top < image.height(); top += options.tile_size) {
    for(std::size_t left = 0; left < image.width(); left += options.tile_size) {
      const std::size_t bottom =
          std::min(top + options.tile_size, image.height());
      const std::size_t right =
          std::min(left + options.tile_size, image.width());
      std::vector<Pixel> pixels;
      pixels.reserve((bottom - top) * (right - left));
      for(std::size_t y = top; y < bottom; ++y) {
        for(std::size_t x = left; x < right; ++x) {
          const astrocfa::RgbPixel pixel = image.pixel(x, y);
          const std::array<double, 3> value = {pixel.r, pixel.g, pixel.b};
          if(std::all_of(value.begin(), value.end(),
                         [](double component) { return std::isfinite(component); })) {
            pixels.push_back(Pixel{
                .luminance =
                    0.2126 * value[0] + 0.7152 * value[1] + 0.0722 * value[2],
                .value = value,
            });
          }
        }
      }
      if(pixels.empty()) {
        continue;
      }
      const std::size_t selected = std::clamp<std::size_t>(
          static_cast<std::size_t>(std::ceil(
              options.sample_quantile * static_cast<double>(pixels.size()))),
          1U, pixels.size());
      std::nth_element(pixels.begin(), pixels.begin() + selected - 1U,
                       pixels.end(), [](const Pixel &left_pixel,
                                        const Pixel &right_pixel) {
                         return left_pixel.luminance < right_pixel.luminance;
                       });
      pixels.resize(selected);
      TileSample sample;
      sample.x = image.width() > 1U
                     ? 2.0 * (0.5 * (left + right - 1U)) /
                               static_cast<double>(image.width() - 1U) -
                           1.0
                     : 0.0;
      sample.y = image.height() > 1U
                     ? 2.0 * (0.5 * (top + bottom - 1U)) /
                               static_cast<double>(image.height() - 1U) -
                           1.0
                     : 0.0;
      for(std::size_t channel = 0; channel < 3; ++channel) {
        std::vector<double> values;
        values.reserve(pixels.size());
        for(const Pixel &pixel : pixels) {
          values.push_back(pixel.value[channel]);
        }
        sample.value[channel] = median(std::move(values));
      }
      samples.push_back(sample);
    }
  }
  return samples;
}

} // namespace

namespace astrocfa {

BackgroundModelResult model_astro_background(
    const RgbImage &linear_camera_rgb, BackgroundModelOptions options) {
  if(options.tile_size < 8U || options.polynomial_degree > 2U ||
     !std::isfinite(options.sample_quantile) ||
     options.sample_quantile <= 0.0 || options.sample_quantile > 0.5 ||
     options.robust_iterations == 0U || options.robust_iterations > 32U ||
     !std::isfinite(options.positive_rejection_sigma) ||
     !std::isfinite(options.negative_rejection_sigma) ||
     options.positive_rejection_sigma <= 0.0 ||
     options.negative_rejection_sigma <= 0.0) {
    throw std::invalid_argument("Invalid astro background model options");
  }
  const std::size_t terms = term_count(options.polynomial_degree);
  std::vector<TileSample> samples = collect_tile_samples(linear_camera_rgb, options);
  while(samples.size() < terms + 2U && options.tile_size > 8U) {
    options.tile_size = std::max<std::size_t>(8U, options.tile_size / 2U);
    samples = collect_tile_samples(linear_camera_rgb, options);
  }
  if(samples.size() < terms + 2U) {
    throw std::invalid_argument(
        "Image has too few background tiles for the selected polynomial");
  }

  std::vector<double> weights(samples.size(), 1.0);
  std::array<std::array<double, maximum_terms>, 3> coefficients{};
  double robust_sigma = 0.0;
  for(std::size_t iteration = 0; iteration < options.robust_iterations;
      ++iteration) {
    for(std::size_t channel = 0; channel < 3; ++channel) {
      coefficients[channel] =
          solve_channel(samples, weights, channel, terms);
    }
    std::vector<double> residuals;
    residuals.reserve(samples.size());
    for(const TileSample &sample : samples) {
      double residual = 0.0;
      for(std::size_t channel = 0; channel < 3; ++channel) {
        residual += sample.value[channel] -
                    evaluate(coefficients[channel], terms, sample.x, sample.y);
      }
      residuals.push_back(residual / 3.0);
    }
    const double center = median(residuals);
    std::vector<double> deviations;
    deviations.reserve(residuals.size());
    for(double residual : residuals) {
      deviations.push_back(std::abs(residual - center));
    }
    robust_sigma = 1.4826 * median(std::move(deviations));
    if(robust_sigma <= 1.0e-12) {
      break;
    }
    for(std::size_t index = 0; index < residuals.size(); ++index) {
      const double centered = residuals[index] - center;
      const double rejection_sigma =
          centered >= 0.0 ? options.positive_rejection_sigma
                          : options.negative_rejection_sigma;
      const double normalized = centered / (rejection_sigma * robust_sigma);
      if(std::abs(normalized) >= 1.0) {
        weights[index] = 1.0e-4;
      } else {
        const double complement = 1.0 - normalized * normalized;
        weights[index] = std::max(1.0e-4, complement * complement);
      }
    }
  }
  for(std::size_t channel = 0; channel < 3; ++channel) {
    coefficients[channel] = solve_channel(samples, weights, channel, terms);
  }

  BackgroundModelResult result{
      .corrected = RgbImage(linear_camera_rgb.width(), linear_camera_rgb.height()),
      .background = RgbImage(linear_camera_rgb.width(), linear_camera_rgb.height()),
  };
  result.stats.effective_tile_size = options.tile_size;
  result.stats.tile_samples = samples.size();
  result.stats.robust_residual_sigma = robust_sigma;
  for(double weight : weights) {
    result.stats.downweighted_samples += weight < 0.25 ? 1U : 0U;
  }

  for(std::size_t channel = 0; channel < 3; ++channel) {
    std::vector<double> modeled_samples;
    modeled_samples.reserve(samples.size());
    for(const TileSample &sample : samples) {
      modeled_samples.push_back(
          evaluate(coefficients[channel], terms, sample.x, sample.y));
    }
    result.stats.preserved_level[channel] = median(modeled_samples);
  }
  if(options.neutralize) {
    const double neutral = median({result.stats.preserved_level[0],
                                   result.stats.preserved_level[1],
                                   result.stats.preserved_level[2]});
    result.stats.preserved_level = {neutral, neutral, neutral};
  }

  std::array<double, 3> minimum = {
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
  };
  std::array<double, 3> maximum = {
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
  };
  for(std::size_t y = 0; y < linear_camera_rgb.height(); ++y) {
    const double normalized_y = linear_camera_rgb.height() > 1U
                                    ? 2.0 * y /
                                              static_cast<double>(
                                                  linear_camera_rgb.height() - 1U) -
                                          1.0
                                    : 0.0;
    for(std::size_t x = 0; x < linear_camera_rgb.width(); ++x) {
      const double normalized_x = linear_camera_rgb.width() > 1U
                                      ? 2.0 * x /
                                                static_cast<double>(
                                                    linear_camera_rgb.width() - 1U) -
                                            1.0
                                      : 0.0;
      const RgbPixel source = linear_camera_rgb.pixel(x, y);
      const std::array<double, 3> source_values = {source.r, source.g, source.b};
      std::array<double, 3> background_values{};
      std::array<double, 3> corrected_values{};
      for(std::size_t channel = 0; channel < 3; ++channel) {
        background_values[channel] = evaluate(
            coefficients[channel], terms, normalized_x, normalized_y);
        corrected_values[channel] =
            source_values[channel] - background_values[channel] +
            result.stats.preserved_level[channel];
        minimum[channel] = std::min(minimum[channel], background_values[channel]);
        maximum[channel] = std::max(maximum[channel], background_values[channel]);
      }
      result.background.set_pixel(
          x, y,
          {.r = static_cast<float>(background_values[0]),
           .g = static_cast<float>(background_values[1]),
           .b = static_cast<float>(background_values[2])});
      result.corrected.set_pixel(
          x, y,
          {.r = static_cast<float>(corrected_values[0]),
           .g = static_cast<float>(corrected_values[1]),
           .b = static_cast<float>(corrected_values[2])});
    }
  }
  for(std::size_t channel = 0; channel < 3; ++channel) {
    result.stats.gradient_peak_to_peak[channel] =
        maximum[channel] - minimum[channel];
  }
  return result;
}

} // namespace astrocfa
