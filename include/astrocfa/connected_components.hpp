#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace astrocfa {

struct BinaryComponentStats {
  std::size_t components = 0;
  std::size_t largest_area = 0;
  std::size_t marked_samples = 0;
};

[[nodiscard]] BinaryComponentStats analyze_binary_components(const std::vector<std::uint8_t> &mask,
                                                             std::size_t width,
                                                             std::size_t height);

} // namespace astrocfa

