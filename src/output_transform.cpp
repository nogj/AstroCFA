#include "astrocfa/output_transform.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

double luminance(astrocfa::RgbPixel pixel) {
  return 0.2126 * pixel.r + 0.7152 * pixel.g + 0.0722 * pixel.b;
}

double percentile(std::vector<double> values, double p) {
  if(values.empty()) {
    return 0.0;
  }
  p = std::clamp(p, 0.0, 100.0);
  const double index = (p / 100.0) * static_cast<double>(values.size() - 1U);
  const auto lower = static_cast<std::size_t>(std::floor(index));
  const auto upper = static_cast<std::size_t>(std::ceil(index));
  std::nth_element(values.begin(), values.begin() + lower, values.end());
  const double lower_value = values[lower];
  if(upper == lower) {
    return lower_value;
  }
  std::nth_element(values.begin(), values.begin() + upper, values.end());
  const double upper_value = values[upper];
  const double fraction = index - static_cast<double>(lower);
  return lower_value * (1.0 - fraction) + upper_value * fraction;
}

struct BaseTransform {
  astrocfa::ToneCurve curve = astrocfa::ToneCurve::linear;
  double distortion = 0.0;
  double local_intensity = 0.0;

  double value(double x) const {
    if(curve == astrocfa::ToneCurve::linear || distortion <= 0.0) {
      return x;
    }
    if(curve == astrocfa::ToneCurve::arcsinh) {
      return std::asinh(distortion * x);
    }
    const double b = local_intensity;
    if(std::abs(b + 1.0) < 1.0e-10) {
      return std::log1p(distortion * x);
    }
    if(b < 0.0) {
      return (1.0 - std::pow(1.0 - b * distortion * x,
                             (b + 1.0) / b)) /
             (distortion * (b + 1.0));
    }
    if(std::abs(b) < 1.0e-10) {
      return 1.0 - std::exp(-distortion * x);
    }
    return 1.0 - std::pow(1.0 + b * distortion * x, -1.0 / b);
  }

  double derivative(double x) const {
    if(curve == astrocfa::ToneCurve::linear || distortion <= 0.0) {
      return 1.0;
    }
    if(curve == astrocfa::ToneCurve::arcsinh) {
      return distortion /
             std::sqrt(1.0 + distortion * distortion * x * x);
    }
    const double b = local_intensity;
    if(std::abs(b + 1.0) < 1.0e-10) {
      return distortion / (1.0 + distortion * x);
    }
    if(b < 0.0) {
      return std::pow(1.0 - b * distortion * x, 1.0 / b);
    }
    if(std::abs(b) < 1.0e-10) {
      return distortion * std::exp(-distortion * x);
    }
    return distortion *
           std::pow(1.0 + b * distortion * x, -(1.0 + b) / b);
  }
};

class FullTransform {
public:
  explicit FullTransform(const astrocfa::AstroToneOptions &options)
      : base_{.curve = options.curve,
              .distortion = std::expm1(options.stretch_factor),
              .local_intensity = options.local_intensity},
        symmetry_(options.symmetry_point), shadows_(options.shadow_protection),
        highlights_(options.highlight_protection) {
    lower_origin_ = lower_segment(0.0);
    range_ = upper_segment(1.0) - lower_origin_;
    if(!std::isfinite(range_) || range_ <= 0.0) {
      throw std::invalid_argument("Tone transform has an invalid output range");
    }
  }

  double operator()(double x) const {
    x = std::clamp(x, 0.0, 1.0);
    double transformed = 0.0;
    if(x < shadows_) {
      transformed = lower_segment(x);
    } else if(x < symmetry_) {
      transformed = -base_.value(symmetry_ - x);
    } else if(x < highlights_) {
      transformed = base_.value(x - symmetry_);
    } else {
      transformed = upper_segment(x);
    }
    return std::clamp((transformed - lower_origin_) / range_, 0.0, 1.0);
  }

private:
  double lower_segment(double x) const {
    const double distance = symmetry_ - shadows_;
    return base_.derivative(distance) * (x - shadows_) -
           base_.value(distance);
  }

  double upper_segment(double x) const {
    const double distance = highlights_ - symmetry_;
    return base_.derivative(distance) * (x - highlights_) +
           base_.value(distance);
  }

  BaseTransform base_;
  double symmetry_ = 0.0;
  double shadows_ = 0.0;
  double highlights_ = 1.0;
  double lower_origin_ = 0.0;
  double range_ = 1.0;
};

} // namespace

namespace astrocfa {

AstroToneResult apply_astro_tone(const RgbImage &linear,
                                 AstroToneOptions options) {
  if(!std::isfinite(options.exposure_ev) || options.exposure_ev < -20.0 ||
     options.exposure_ev > 20.0 || !std::isfinite(options.black_percentile) ||
     !std::isfinite(options.white_percentile) ||
     options.black_percentile < 0.0 || options.white_percentile > 100.0 ||
     options.black_percentile >= options.white_percentile ||
     !std::isfinite(options.black_point) ||
     !std::isfinite(options.white_point) ||
     !std::isfinite(options.stretch_factor) || options.stretch_factor < 0.0 ||
     options.stretch_factor > 20.0 ||
     !std::isfinite(options.local_intensity) ||
     options.local_intensity < -5.0 || options.local_intensity > 15.0 ||
     !std::isfinite(options.symmetry_point) ||
     !std::isfinite(options.shadow_protection) ||
     !std::isfinite(options.highlight_protection) ||
     options.shadow_protection < 0.0 ||
     options.shadow_protection > options.symmetry_point ||
     options.symmetry_point > options.highlight_protection ||
     options.highlight_protection > 1.0 ||
     !std::isfinite(options.saturation) || options.saturation < 0.0 ||
     options.saturation > 4.0) {
    throw std::invalid_argument("Invalid astro tone options");
  }

  const double exposure = std::exp2(options.exposure_ev);
  std::vector<double> luminance_samples;
  luminance_samples.reserve(linear.width() * linear.height());
  double maximum_input_luminance = 0.0;
  for(std::size_t y = 0; y < linear.height(); ++y) {
    for(std::size_t x = 0; x < linear.width(); ++x) {
      const double value = luminance(linear.pixel(x, y));
      if(std::isfinite(value)) {
        luminance_samples.push_back(value);
        maximum_input_luminance =
            std::max(maximum_input_luminance, exposure * value);
      }
    }
  }
  double black = options.black_point;
  double white = options.white_point;
  if(options.auto_levels) {
    black = percentile(luminance_samples, options.black_percentile);
    white = percentile(luminance_samples, options.white_percentile);
  }
  if(!std::isfinite(black) || !std::isfinite(white) || white <= black) {
    throw std::invalid_argument("Tone white point must be greater than black point");
  }

  const FullTransform transform(options);
  AstroToneResult result{
      .image = RgbImage(linear.width(), linear.height()),
      .stats = AstroToneStats{
          .black_point = black,
          .white_point = white,
          .maximum_input_luminance = maximum_input_luminance,
      },
  };
  const double dynamic_range = white - black;
  for(std::size_t y = 0; y < linear.height(); ++y) {
    for(std::size_t x = 0; x < linear.width(); ++x) {
      const RgbPixel source = linear.pixel(x, y);
      std::array<double, 3> normalized = {
          std::max(0.0, (exposure * source.r - black) / dynamic_range),
          std::max(0.0, (exposure * source.g - black) / dynamic_range),
          std::max(0.0, (exposure * source.b - black) / dynamic_range),
      };
      const double source_luminance =
          0.2126 * normalized[0] + 0.7152 * normalized[1] +
          0.0722 * normalized[2];
      result.stats.shadow_clipped_pixels += source_luminance <= 0.0 ? 1U : 0U;
      result.stats.highlight_clipped_pixels += source_luminance >= 1.0 ? 1U : 0U;
      const double mapped_luminance = transform(source_luminance);
      const double scale = source_luminance > 1.0e-15
                               ? mapped_luminance / source_luminance
                               : 0.0;
      std::array<double, 3> output{};
      for(std::size_t channel = 0; channel < 3; ++channel) {
        const double color_preserved = normalized[channel] * scale;
        output[channel] = mapped_luminance +
                          options.saturation *
                              (color_preserved - mapped_luminance);
      }
      const double maximum = std::max({output[0], output[1], output[2]});
      const double minimum = std::min({output[0], output[1], output[2]});
      if(maximum > 1.0 || minimum < 0.0) {
        const double chroma_scale =
            maximum > 1.0
                ? std::min(1.0, (1.0 - mapped_luminance) /
                                    std::max(1.0e-15, maximum - mapped_luminance))
                : std::min(1.0, mapped_luminance /
                                    std::max(1.0e-15, mapped_luminance - minimum));
        for(double &component : output) {
          component = mapped_luminance +
                      chroma_scale * (component - mapped_luminance);
        }
        result.stats.gamut_compressed_pixels += 1U;
      }
      for(double &component : output) {
        component = std::clamp(component, 0.0, 1.0);
        result.stats.maximum_output_component =
            std::max(result.stats.maximum_output_component, component);
      }
      result.image.set_pixel(
          x, y,
          {.r = static_cast<float>(output[0]),
           .g = static_cast<float>(output[1]),
           .b = static_cast<float>(output[2])});
    }
  }
  return result;
}

RgbImage make_astro_preview(const RgbImage &linear, AstroPreviewOptions options) {
  if(!std::isfinite(options.arcsinh_strength) ||
     options.arcsinh_strength <= 0.0) {
    throw std::invalid_argument("Astro preview arcsinh strength must be positive");
  }
  return apply_astro_tone(
             linear, AstroToneOptions{
                         .curve = ToneCurve::arcsinh,
                         .black_percentile = options.black_percentile,
                         .white_percentile = options.white_percentile,
                         .stretch_factor = std::log1p(options.arcsinh_strength),
                     })
      .image;
}

} // namespace astrocfa
