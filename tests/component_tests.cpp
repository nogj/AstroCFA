#include "astrocfa/connected_components.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char *message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void detects_8_connected_components() {
  std::vector<std::uint8_t> mask = {
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
  };

  const astrocfa::BinaryComponentStats stats =
      astrocfa::analyze_binary_components(mask, 4, 4);
  require(stats.components == 1, "Diagonal samples should be connected");
  require(stats.largest_area == 4, "Largest component should include diagonal chain");
  require(stats.marked_samples == 4, "Marked sample count should match mask");
}

void separates_isolated_samples() {
  std::vector<std::uint8_t> mask = {
      1, 0, 0, 1,
      0, 0, 0, 0,
      0, 0, 0, 0,
      1, 0, 0, 1,
  };

  const astrocfa::BinaryComponentStats stats =
      astrocfa::analyze_binary_components(mask, 4, 4);
  require(stats.components == 4, "Corner samples should be separate components");
  require(stats.largest_area == 1, "Each isolated sample has area one");
  require(stats.marked_samples == 4, "Marked sample count should match mask");
}

} // namespace

int main() {
  try {
    detects_8_connected_components();
    separates_isolated_samples();
  } catch(const std::exception &error) {
    std::cerr << "component_tests failed: " << error.what() << "\n";
    return 1;
  }

  return 0;
}

