#include "astrocfa/cfa.hpp"

#include <algorithm>
#include <cmath>

namespace astrocfa {

CfaFrame::CfaFrame(std::size_t width, std::size_t height, BayerPattern pattern)
    : width_(width), height_(height), pattern_(pattern), samples_(width * height) {}

float CfaFrame::sample(std::size_t x, std::size_t y) const {
  return samples_[offset(x, y)].value;
}

CfaSample CfaFrame::sample_info(std::size_t x, std::size_t y) const {
  return samples_[offset(x, y)];
}

void CfaFrame::set_sample(std::size_t x, std::size_t y, float value) {
  samples_[offset(x, y)].value = value;
}

void CfaFrame::set_sample(std::size_t x, std::size_t y, CfaSample sample) {
  samples_[offset(x, y)] = sample;
}

void CfaFrame::set_valid(std::size_t x, std::size_t y, bool valid) {
  samples_[offset(x, y)].valid = valid;
}

void CfaFrame::set_clipped(std::size_t x, std::size_t y, bool clipped) {
  samples_[offset(x, y)].clipped = clipped;
}

std::size_t CfaFrame::offset(std::size_t x, std::size_t y) const {
  if(x >= width_ || y >= height_) {
    throw std::out_of_range("CFA coordinate out of range");
  }
  return y * width_ + x;
}

RgbImage::RgbImage(std::size_t width, std::size_t height)
    : width_(width), height_(height), pixels_(width * height) {}

RgbPixel RgbImage::pixel(std::size_t x, std::size_t y) const {
  return pixels_[offset(x, y)];
}

void RgbImage::set_pixel(std::size_t x, std::size_t y, RgbPixel value) {
  pixels_[offset(x, y)] = value;
}

std::size_t RgbImage::offset(std::size_t x, std::size_t y) const {
  if(x >= width_ || y >= height_) {
    throw std::out_of_range("RGB coordinate out of range");
  }
  return y * width_ + x;
}

float remosaic_pixel(RgbPixel pixel, CfaColor color) {
  switch(color) {
  case CfaColor::red:
    return pixel.r;
  case CfaColor::green1:
  case CfaColor::green2:
    return pixel.g;
  case CfaColor::blue:
    return pixel.b;
  }

  return 0.0F;
}

CfaFrame remosaic(const RgbImage &image, BayerPattern pattern) {
  CfaFrame cfa(image.width(), image.height(), pattern);
  for(std::size_t y = 0; y < image.height(); ++y) {
    for(std::size_t x = 0; x < image.width(); ++x) {
      cfa.set_sample(x, y, remosaic_pixel(image.pixel(x, y), pattern.at(x, y)));
    }
  }
  return cfa;
}

RemosaicResidual compute_remosaic_residual(const CfaFrame &measured,
                                           const RgbImage &reconstructed) {
  if(measured.width() != reconstructed.width() || measured.height() != reconstructed.height()) {
    throw std::invalid_argument("Measured CFA and reconstructed RGB dimensions differ");
  }

  double absolute_sum = 0.0;
  double square_sum = 0.0;
  float maximum_absolute = 0.0F;
  std::size_t samples = 0;

  const BayerPattern &pattern = measured.pattern();
  for(std::size_t y = 0; y < measured.height(); ++y) {
    for(std::size_t x = 0; x < measured.width(); ++x) {
      const float predicted = remosaic_pixel(reconstructed.pixel(x, y), pattern.at(x, y));
      const CfaSample measured_sample = measured.sample_info(x, y);
      if(!measured_sample.valid || measured_sample.clipped) {
        continue;
      }

      const float residual = predicted - measured_sample.value;
      const float absolute = std::abs(residual);
      absolute_sum += absolute;
      square_sum += static_cast<double>(residual) * static_cast<double>(residual);
      maximum_absolute = std::max(maximum_absolute, absolute);
      ++samples;
    }
  }

  RemosaicResidual residual;
  residual.samples = samples;
  if(samples > 0) {
    residual.mean_absolute = absolute_sum / static_cast<double>(samples);
    residual.root_mean_square = std::sqrt(square_sum / static_cast<double>(samples));
    residual.maximum_absolute = maximum_absolute;
  }
  return residual;
}

NoiseWeightedResidual
compute_noise_weighted_remosaic_residual(const CfaFrame &measured,
                                         const RgbImage &reconstructed,
                                         const NoiseModel &noise_model) {
  if(measured.width() != reconstructed.width() || measured.height() != reconstructed.height()) {
    throw std::invalid_argument("Measured CFA and reconstructed RGB dimensions differ");
  }

  double chi_square = 0.0;
  double normalized_absolute_sum = 0.0;
  double max_normalized_absolute = 0.0;
  std::size_t samples = 0;

  const BayerPattern &pattern = measured.pattern();
  for(std::size_t y = 0; y < measured.height(); ++y) {
    for(std::size_t x = 0; x < measured.width(); ++x) {
      const CfaSample measured_sample = measured.sample_info(x, y);
      if(!measured_sample.valid || measured_sample.clipped) {
        continue;
      }

      const float predicted = remosaic_pixel(reconstructed.pixel(x, y), pattern.at(x, y));
      const double residual = static_cast<double>(predicted) - measured_sample.value;
      const NoiseEstimate noise = estimate_noise(measured_sample.value, noise_model);
      const double normalized = residual / noise.sigma;
      const double absolute = std::abs(normalized);
      chi_square += normalized * normalized;
      normalized_absolute_sum += absolute;
      max_normalized_absolute = std::max(max_normalized_absolute, absolute);
      ++samples;
    }
  }

  NoiseWeightedResidual residual;
  residual.samples = samples;
  residual.chi_square = chi_square;
  if(samples > 0) {
    residual.reduced_chi_square = chi_square / static_cast<double>(samples);
    residual.mean_normalized_absolute =
        normalized_absolute_sum / static_cast<double>(samples);
    residual.max_normalized_absolute = max_normalized_absolute;
  }
  return residual;
}

} // namespace astrocfa
