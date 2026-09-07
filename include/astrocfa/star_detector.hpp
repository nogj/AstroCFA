#pragma once

#include "astrocfa/cfa.hpp"

#include <cstddef>

namespace astrocfa {

struct StarDetectionStats {
  std::size_t candidates = 0;
  std::size_t largest_area = 0;
  double background_mean = 0.0;
  double background_sigma = 0.0;
  double threshold = 0.0;
  float brightest = 0.0F;
};

[[nodiscard]] StarDetectionStats detect_star_candidates(const CfaFrame &cfa,
                                                        double sigma_threshold = 5.0);

} // namespace astrocfa

