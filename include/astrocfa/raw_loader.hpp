#pragma once

#include "astrocfa/cfa.hpp"
#include "astrocfa/raw_inspector.hpp"

#include <string>

namespace astrocfa {

struct LinearRawFrame {
  RawInspection inspection;
  CfaFrame cfa;
};

[[nodiscard]] LinearRawFrame load_linear_cfa_file(const std::string &path);

} // namespace astrocfa

