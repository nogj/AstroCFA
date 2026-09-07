#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/drizzle.hpp"

#include <cstddef>

namespace astrocfa {

struct RegistrationOptions {
  int max_shift_pixels = 64;
  double sigma_threshold = 4.0;
  std::size_t sample_stride = 16;
};

struct RegistrationResult {
  SubpixelOffset offset;
  double score = 0.0;
  int proxy_dx = 0;
  int proxy_dy = 0;
  std::size_t matched_samples = 0;
};

[[nodiscard]] RegistrationResult estimate_integer_star_translation(
    const CfaFrame &reference, const CfaFrame &moving,
    RegistrationOptions options = {});

} // namespace astrocfa
