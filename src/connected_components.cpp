#include "astrocfa/connected_components.hpp"

#include <algorithm>
#include <stdexcept>

namespace astrocfa {

BinaryComponentStats analyze_binary_components(const std::vector<std::uint8_t> &mask,
                                               std::size_t width,
                                               std::size_t height) {
  if(mask.size() != width * height) {
    throw std::invalid_argument("Binary mask dimensions do not match mask size");
  }

  BinaryComponentStats stats;
  std::vector<std::uint8_t> visited(mask.size(), 0);
  std::vector<std::size_t> stack;

  const auto offset = [width](std::size_t x, std::size_t y) { return y * width + x; };

  for(std::size_t y = 0; y < height; ++y) {
    for(std::size_t x = 0; x < width; ++x) {
      const std::size_t start = offset(x, y);
      if(mask[start] == 0 || visited[start] != 0) {
        continue;
      }

      std::size_t area = 0;
      stats.components += 1;
      stack.clear();
      stack.push_back(start);
      visited[start] = 1;

      while(!stack.empty()) {
        const std::size_t current = stack.back();
        stack.pop_back();
        area += 1;

        const std::size_t cy = current / width;
        const std::size_t cx = current - cy * width;
        const std::size_t x0 = cx == 0 ? 0 : cx - 1;
        const std::size_t y0 = cy == 0 ? 0 : cy - 1;
        const std::size_t x1 = std::min(width - 1, cx + 1);
        const std::size_t y1 = std::min(height - 1, cy + 1);

        for(std::size_t ny = y0; ny <= y1; ++ny) {
          for(std::size_t nx = x0; nx <= x1; ++nx) {
            const std::size_t neighbor = offset(nx, ny);
            if(mask[neighbor] == 0 || visited[neighbor] != 0) {
              continue;
            }
            visited[neighbor] = 1;
            stack.push_back(neighbor);
          }
        }
      }

      stats.marked_samples += area;
      stats.largest_area = std::max(stats.largest_area, area);
    }
  }

  return stats;
}

} // namespace astrocfa

