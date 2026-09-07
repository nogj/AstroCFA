#include "astrocfa/demosaic.hpp"

#include <algorithm>
#include <cmath>

namespace {

int channel_for_phase(astrocfa::CfaColor color) {
  switch(color) {
  case astrocfa::CfaColor::red:
    return 0;
  case astrocfa::CfaColor::green1:
  case astrocfa::CfaColor::green2:
    return 1;
  case astrocfa::CfaColor::blue:
    return 2;
  }
  return 1;
}

float measured_channel_value(const astrocfa::CfaSample &sample) {
  return sample.valid ? sample.value : 0.0F;
}

float interpolate_channel(const astrocfa::CfaFrame &cfa, std::size_t x, std::size_t y,
                          int target_channel) {
  double sum = 0.0;
  std::size_t count = 0;

  const auto width = static_cast<int>(cfa.width());
  const auto height = static_cast<int>(cfa.height());
  const int ix = static_cast<int>(x);
  const int iy = static_cast<int>(y);

  for(int dy = -1; dy <= 1; ++dy) {
    for(int dx = -1; dx <= 1; ++dx) {
      const int nx = ix + dx;
      const int ny = iy + dy;
      if(nx < 0 || ny < 0 || nx >= width || ny >= height) {
        continue;
      }

      const auto phase = cfa.pattern().at(static_cast<std::size_t>(nx),
                                          static_cast<std::size_t>(ny));
      if(channel_for_phase(phase) != target_channel) {
        continue;
      }

      const astrocfa::CfaSample sample =
          cfa.sample_info(static_cast<std::size_t>(nx), static_cast<std::size_t>(ny));
      if(!sample.valid) {
        continue;
      }

      sum += measured_channel_value(sample);
      count += 1;
    }
  }

  if(count > 0) {
    return static_cast<float>(sum / static_cast<double>(count));
  }

  return measured_channel_value(cfa.sample_info(x, y));
}

float interpolate_residual(const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &green_image,
                           std::size_t x, std::size_t y, int target_channel) {
  double sum = 0.0;
  double weight_sum = 0.0;

  const auto width = static_cast<int>(cfa.width());
  const auto height = static_cast<int>(cfa.height());
  const int ix = static_cast<int>(x);
  const int iy = static_cast<int>(y);

  for(int dy = -2; dy <= 2; ++dy) {
    for(int dx = -2; dx <= 2; ++dx) {
      const int nx = ix + dx;
      const int ny = iy + dy;
      if(nx < 0 || ny < 0 || nx >= width || ny >= height) {
        continue;
      }

      const auto ux = static_cast<std::size_t>(nx);
      const auto uy = static_cast<std::size_t>(ny);
      if(channel_for_phase(cfa.pattern().at(ux, uy)) != target_channel) {
        continue;
      }

      const astrocfa::CfaSample sample = cfa.sample_info(ux, uy);
      if(!sample.valid || sample.clipped) {
        continue;
      }

      const float green = green_image.pixel(ux, uy).g;
      const double distance2 = static_cast<double>(dx * dx + dy * dy);
      const double weight = 1.0 / (1.0 + distance2);
      sum += (sample.value - green) * weight;
      weight_sum += weight;
    }
  }

  return weight_sum > 0.0 ? static_cast<float>(sum / weight_sum) : 0.0F;
}

float safe_sample(const astrocfa::CfaFrame &cfa, int x, int y) {
  const int clamped_x = std::clamp(x, 0, static_cast<int>(cfa.width()) - 1);
  const int clamped_y = std::clamp(y, 0, static_cast<int>(cfa.height()) - 1);
  return measured_channel_value(
      cfa.sample_info(static_cast<std::size_t>(clamped_x), static_cast<std::size_t>(clamped_y)));
}

float clamp_unit(float value) {
  return std::clamp(value, 0.0F, 1.25F);
}

float channel_value(astrocfa::RgbPixel pixel, int channel) {
  if(channel == 0) {
    return pixel.r;
  }
  if(channel == 2) {
    return pixel.b;
  }
  return pixel.g;
}

void set_channel_value(astrocfa::RgbPixel &pixel, int channel, float value) {
  if(channel == 0) {
    pixel.r = value;
  } else if(channel == 2) {
    pixel.b = value;
  } else {
    pixel.g = value;
  }
}

double luminance_gradient(const astrocfa::RgbImage &image, std::size_t x0, std::size_t y0,
                          std::size_t x1, std::size_t y1) {
  const auto luma = [](astrocfa::RgbPixel pixel) {
    return 0.2126 * pixel.r + 0.7152 * pixel.g + 0.0722 * pixel.b;
  };
  return std::abs(luma(image.pixel(x0, y0)) - luma(image.pixel(x1, y1)));
}

float refined_chroma_at(const astrocfa::CfaFrame &cfa, const astrocfa::RgbImage &current,
                        const astrocfa::FrequencyRiskMap &risk_map, std::size_t x,
                        std::size_t y, int chroma_channel,
                        const astrocfa::InverseRefinementOptions &options) {
  const int measured_channel = channel_for_phase(cfa.pattern().at(x, y));
  const astrocfa::RgbPixel center = current.pixel(x, y);
  const double center_chroma = channel_value(center, chroma_channel) - center.g;

  if(measured_channel == chroma_channel) {
    return static_cast<float>(center_chroma);
  }

  const auto width = static_cast<int>(cfa.width());
  const auto height = static_cast<int>(cfa.height());
  const int ix = static_cast<int>(x);
  const int iy = static_cast<int>(y);
  double weighted_sum = 0.0;
  double weight_sum = 0.0;

  const int offsets[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  for(const auto &offset : offsets) {
    const int nx = ix + offset[0];
    const int ny = iy + offset[1];
    if(nx < 0 || ny < 0 || nx >= width || ny >= height) {
      continue;
    }

    const auto ux = static_cast<std::size_t>(nx);
    const auto uy = static_cast<std::size_t>(ny);
    const astrocfa::RgbPixel neighbor = current.pixel(ux, uy);
    const double neighbor_chroma = channel_value(neighbor, chroma_channel) - neighbor.g;
    const double edge =
        luminance_gradient(current, x, y, ux, uy) * options.edge_sensitivity;
    const double edge_weight = 1.0 / (1.0 + edge * edge);
    const int neighbor_measured = channel_for_phase(cfa.pattern().at(ux, uy));
    const double anchor_weight = neighbor_measured == chroma_channel ? 2.0 : 1.0;
    const astrocfa::CfaSample sample = cfa.sample_info(ux, uy);
    const double sample_weight = sample.valid && !sample.clipped ? 1.0 : 0.25;
    const double weight = edge_weight * anchor_weight * sample_weight;

    weighted_sum += neighbor_chroma * weight;
    weight_sum += weight;
  }

  if(weight_sum <= 0.0) {
    return static_cast<float>(center_chroma);
  }

  const double neighbor_average = weighted_sum / weight_sum;
  const double alias = std::clamp(risk_map.pixel_risk(x, y), 0.0, 1.0);
  const double smoothness =
      std::clamp(options.chroma_smoothness + alias * options.alias_suppression, 0.0, 1.0);
  return static_cast<float>((1.0 - smoothness) * center_chroma +
                            smoothness * neighbor_average);
}

void restore_measured_channel(const astrocfa::CfaFrame &cfa, std::size_t x, std::size_t y,
                              astrocfa::RgbPixel &pixel) {
  const int measured_channel = channel_for_phase(cfa.pattern().at(x, y));
  const astrocfa::CfaSample measured = cfa.sample_info(x, y);
  set_channel_value(pixel, measured_channel, measured_channel_value(measured));
}

} // namespace

namespace astrocfa {

RgbImage demosaic_bilinear_baseline(const CfaFrame &cfa) {
  RgbImage rgb(cfa.width(), cfa.height());

  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const CfaColor phase = cfa.pattern().at(x, y);
      const int measured_channel = channel_for_phase(phase);
      const CfaSample measured = cfa.sample_info(x, y);

      RgbPixel pixel;
      pixel.r = measured_channel == 0 ? measured_channel_value(measured)
                                      : interpolate_channel(cfa, x, y, 0);
      pixel.g = measured_channel == 1 ? measured_channel_value(measured)
                                      : interpolate_channel(cfa, x, y, 1);
      pixel.b = measured_channel == 2 ? measured_channel_value(measured)
                                      : interpolate_channel(cfa, x, y, 2);
      rgb.set_pixel(x, y, pixel);
    }
  }

  return rgb;
}

RgbImage demosaic_malvar_baseline(const CfaFrame &cfa) {
  RgbImage rgb = demosaic_bilinear_baseline(cfa);

  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const int ix = static_cast<int>(x);
      const int iy = static_cast<int>(y);
      const CfaColor phase = cfa.pattern().at(x, y);
      const int measured_channel = channel_for_phase(phase);
      const float center = safe_sample(cfa, ix, iy);
      RgbPixel pixel = rgb.pixel(x, y);

      if(measured_channel == 0 || measured_channel == 2) {
        const float axial_average =
            0.25F * (safe_sample(cfa, ix - 1, iy) + safe_sample(cfa, ix + 1, iy) +
                     safe_sample(cfa, ix, iy - 1) + safe_sample(cfa, ix, iy + 1));
        const float laplacian =
            center - 0.25F * (safe_sample(cfa, ix - 2, iy) + safe_sample(cfa, ix + 2, iy) +
                              safe_sample(cfa, ix, iy - 2) + safe_sample(cfa, ix, iy + 2));
        pixel.g = clamp_unit(axial_average + 0.5F * laplacian);
      } else {
        const bool horizontal_red =
            channel_for_phase(cfa.pattern().at(x > 0 ? x - 1 : x, y)) == 0 ||
            channel_for_phase(cfa.pattern().at(std::min(cfa.width() - 1, x + 1), y)) == 0;
        const int primary_channel = horizontal_red ? 0 : 2;
        const int secondary_channel = horizontal_red ? 2 : 0;
        const float horizontal_average =
            0.5F * (safe_sample(cfa, ix - 1, iy) + safe_sample(cfa, ix + 1, iy));
        const float vertical_average =
            0.5F * (safe_sample(cfa, ix, iy - 1) + safe_sample(cfa, ix, iy + 1));
        const float horizontal_laplacian =
            center - 0.5F * (safe_sample(cfa, ix - 2, iy) + safe_sample(cfa, ix + 2, iy));
        const float vertical_laplacian =
            center - 0.5F * (safe_sample(cfa, ix, iy - 2) + safe_sample(cfa, ix, iy + 2));

        if(primary_channel == 0) {
          pixel.r = clamp_unit(horizontal_average + 0.5F * horizontal_laplacian);
          pixel.b = clamp_unit(vertical_average + 0.5F * vertical_laplacian);
        } else {
          pixel.b = clamp_unit(horizontal_average + 0.5F * horizontal_laplacian);
          pixel.r = clamp_unit(vertical_average + 0.5F * vertical_laplacian);
        }

        (void)secondary_channel;
      }

      if(measured_channel == 0) {
        const float diagonal_average =
            0.25F * (safe_sample(cfa, ix - 1, iy - 1) + safe_sample(cfa, ix + 1, iy - 1) +
                     safe_sample(cfa, ix - 1, iy + 1) + safe_sample(cfa, ix + 1, iy + 1));
        const float diagonal_laplacian =
            center - 0.25F * (safe_sample(cfa, ix - 2, iy) + safe_sample(cfa, ix + 2, iy) +
                              safe_sample(cfa, ix, iy - 2) + safe_sample(cfa, ix, iy + 2));
        pixel.b = clamp_unit(diagonal_average + 0.75F * diagonal_laplacian);
      } else if(measured_channel == 2) {
        const float diagonal_average =
            0.25F * (safe_sample(cfa, ix - 1, iy - 1) + safe_sample(cfa, ix + 1, iy - 1) +
                     safe_sample(cfa, ix - 1, iy + 1) + safe_sample(cfa, ix + 1, iy + 1));
        const float diagonal_laplacian =
            center - 0.25F * (safe_sample(cfa, ix - 2, iy) + safe_sample(cfa, ix + 2, iy) +
                              safe_sample(cfa, ix, iy - 2) + safe_sample(cfa, ix, iy + 2));
        pixel.r = clamp_unit(diagonal_average + 0.75F * diagonal_laplacian);
      }

      const CfaSample measured = cfa.sample_info(x, y);
      if(measured_channel == 0) {
        pixel.r = measured_channel_value(measured);
      } else if(measured_channel == 1) {
        pixel.g = measured_channel_value(measured);
      } else {
        pixel.b = measured_channel_value(measured);
      }

      rgb.set_pixel(x, y, pixel);
    }
  }

  return rgb;
}

RgbImage demosaic_frequency_guided(const CfaFrame &cfa,
                                   FrequencyGuidedDemosaicOptions options) {
  RgbImage rgb = demosaic_residual_interpolation(cfa);
  const FrequencyRiskMap risk_map = build_frequency_risk_map(cfa, options.frequency);

  for(std::size_t y = 0; y < rgb.height(); ++y) {
    for(std::size_t x = 0; x < rgb.width(); ++x) {
      const CfaColor phase = cfa.pattern().at(x, y);
      const int measured_channel = channel_for_phase(phase);
      RgbPixel pixel = rgb.pixel(x, y);
      const float luma_anchor = pixel.g;
      const double suppression =
          std::clamp(risk_map.pixel_risk(x, y) * options.max_chroma_suppression, 0.0, 1.0);

      if(measured_channel != 0) {
        pixel.r =
            static_cast<float>((1.0 - suppression) * pixel.r + suppression * luma_anchor);
      }
      if(measured_channel != 2) {
        pixel.b =
            static_cast<float>((1.0 - suppression) * pixel.b + suppression * luma_anchor);
      }
      rgb.set_pixel(x, y, pixel);
    }
  }

  return rgb;
}

RgbImage demosaic_residual_interpolation(const CfaFrame &cfa) {
  RgbImage rgb = demosaic_malvar_baseline(cfa);

  for(std::size_t y = 0; y < cfa.height(); ++y) {
    for(std::size_t x = 0; x < cfa.width(); ++x) {
      const int measured_channel = channel_for_phase(cfa.pattern().at(x, y));
      const CfaSample measured = cfa.sample_info(x, y);
      RgbPixel pixel = rgb.pixel(x, y);
      const float green = pixel.g;

      if(measured_channel == 0) {
        pixel.r = measured_channel_value(measured);
      } else {
        pixel.r = clamp_unit(green + interpolate_residual(cfa, rgb, x, y, 0));
      }

      if(measured_channel == 2) {
        pixel.b = measured_channel_value(measured);
      } else {
        pixel.b = clamp_unit(green + interpolate_residual(cfa, rgb, x, y, 2));
      }

      if(measured_channel == 1) {
        pixel.g = measured_channel_value(measured);
      }

      rgb.set_pixel(x, y, pixel);
    }
  }

  return rgb;
}

RgbImage demosaic_inverse_refine(const CfaFrame &cfa, InverseRefinementOptions options) {
  RgbImage current = demosaic_frequency_guided(
      cfa, FrequencyGuidedDemosaicOptions{.frequency = options.frequency});
  const FrequencyRiskMap risk_map = build_frequency_risk_map(cfa, options.frequency);
  const int iterations = std::max(0, options.iterations);

  for(int iteration = 0; iteration < iterations; ++iteration) {
    RgbImage next = current;
    for(std::size_t y = 0; y < cfa.height(); ++y) {
      for(std::size_t x = 0; x < cfa.width(); ++x) {
        RgbPixel pixel = current.pixel(x, y);
        const int measured_channel = channel_for_phase(cfa.pattern().at(x, y));
        if(measured_channel != 0) {
          const float rg =
              refined_chroma_at(cfa, current, risk_map, x, y, 0, options);
          pixel.r = clamp_unit(pixel.g + rg);
        }
        if(measured_channel != 2) {
          const float bg =
              refined_chroma_at(cfa, current, risk_map, x, y, 2, options);
          pixel.b = clamp_unit(pixel.g + bg);
        }
        restore_measured_channel(cfa, x, y, pixel);
        next.set_pixel(x, y, pixel);
      }
    }
    current = next;
  }

  return current;
}

DemosaicResult reconstruct_baseline(const CfaFrame &cfa, const NoiseModel &noise_model) {
  RgbImage rgb = demosaic_bilinear_baseline(cfa);
  return DemosaicResult{
      .image = rgb,
      .residual = compute_remosaic_residual(cfa, rgb),
      .noise_weighted_residual =
          compute_noise_weighted_remosaic_residual(cfa, rgb, noise_model),
  };
}

DemosaicResult reconstruct_malvar_baseline(const CfaFrame &cfa,
                                           const NoiseModel &noise_model) {
  RgbImage rgb = demosaic_malvar_baseline(cfa);
  return DemosaicResult{
      .image = rgb,
      .residual = compute_remosaic_residual(cfa, rgb),
      .noise_weighted_residual =
          compute_noise_weighted_remosaic_residual(cfa, rgb, noise_model),
  };
}

DemosaicResult reconstruct_residual_interpolation(const CfaFrame &cfa,
                                                  const NoiseModel &noise_model) {
  RgbImage rgb = demosaic_residual_interpolation(cfa);
  return DemosaicResult{
      .image = rgb,
      .residual = compute_remosaic_residual(cfa, rgb),
      .noise_weighted_residual =
          compute_noise_weighted_remosaic_residual(cfa, rgb, noise_model),
  };
}

DemosaicResult reconstruct_frequency_guided(const CfaFrame &cfa,
                                            const NoiseModel &noise_model,
                                            FrequencyGuidedDemosaicOptions options) {
  RgbImage rgb = demosaic_frequency_guided(cfa, options);
  return DemosaicResult{
      .image = rgb,
      .residual = compute_remosaic_residual(cfa, rgb),
      .noise_weighted_residual =
          compute_noise_weighted_remosaic_residual(cfa, rgb, noise_model),
  };
}

DemosaicResult reconstruct_inverse_refine(const CfaFrame &cfa,
                                          const NoiseModel &noise_model,
                                          InverseRefinementOptions options) {
  RgbImage rgb = demosaic_inverse_refine(cfa, options);
  return DemosaicResult{
      .image = rgb,
      .residual = compute_remosaic_residual(cfa, rgb),
      .noise_weighted_residual =
          compute_noise_weighted_remosaic_residual(cfa, rgb, noise_model),
  };
}

DemosaicQuality analyze_demosaic_quality(const RgbImage &image) {
  DemosaicQuality quality;
  if(image.width() < 3 || image.height() < 3) {
    return quality;
  }

  for(std::size_t y = 1; y + 1 < image.height(); ++y) {
    for(std::size_t x = 1; x + 1 < image.width(); ++x) {
      const RgbPixel center = image.pixel(x, y);
      const double center_rg = static_cast<double>(center.r) - center.g;
      const double center_bg = static_cast<double>(center.b) - center.g;
      double roughness = 0.0;

      const std::size_t xs[4] = {x - 1, x + 1, x, x};
      const std::size_t ys[4] = {y, y, y - 1, y + 1};
      for(std::size_t i = 0; i < 4; ++i) {
        const RgbPixel neighbor = image.pixel(xs[i], ys[i]);
        roughness +=
            std::abs(center_rg - (static_cast<double>(neighbor.r) - neighbor.g));
        roughness +=
            std::abs(center_bg - (static_cast<double>(neighbor.b) - neighbor.g));
      }
      roughness *= 0.25;

      quality.mean_chroma_roughness += roughness;
      quality.max_chroma_roughness = std::max(quality.max_chroma_roughness, roughness);
      quality.samples += 1;
    }
  }

  if(quality.samples > 0) {
    quality.mean_chroma_roughness /= static_cast<double>(quality.samples);
  }
  return quality;
}

DemosaicQuality analyze_demosaic_quality(const RgbImage &image, const CfaFrame &cfa) {
  if(image.width() != cfa.width() || image.height() != cfa.height()) {
    throw std::invalid_argument("Demosaic quality image and CFA dimensions differ");
  }

  DemosaicQuality quality = analyze_demosaic_quality(image);
  double interpolated_chroma_sum = 0.0;

  for(std::size_t y = 0; y < image.height(); ++y) {
    for(std::size_t x = 0; x < image.width(); ++x) {
      const int measured_channel = channel_for_phase(cfa.pattern().at(x, y));
      const RgbPixel pixel = image.pixel(x, y);

      if(measured_channel != 0) {
        const double chroma = std::abs(static_cast<double>(pixel.r) - pixel.g);
        interpolated_chroma_sum += chroma;
        quality.max_interpolated_chroma =
            std::max(quality.max_interpolated_chroma, chroma);
        quality.interpolated_chroma_samples += 1;
      }
      if(measured_channel != 2) {
        const double chroma = std::abs(static_cast<double>(pixel.b) - pixel.g);
        interpolated_chroma_sum += chroma;
        quality.max_interpolated_chroma =
            std::max(quality.max_interpolated_chroma, chroma);
        quality.interpolated_chroma_samples += 1;
      }
    }
  }

  if(quality.interpolated_chroma_samples > 0) {
    quality.mean_interpolated_chroma =
        interpolated_chroma_sum / static_cast<double>(quality.interpolated_chroma_samples);
  }

  return quality;
}

} // namespace astrocfa
