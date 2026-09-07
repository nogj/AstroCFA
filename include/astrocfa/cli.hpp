#pragma once

#include "astrocfa/version.hpp"

#include <ostream>
#include <string_view>

namespace astrocfa {

inline void print_cli_help(std::ostream &out, std::string_view executable) {
  out << "AstroCFA " << version << "\n"
      << tagline << "\n\n"
      << "Usage:\n"
      << "  " << executable << " --help\n"
      << "  " << executable << " --version\n"
      << "  " << executable
      << " inspect <input-raw-or-dng> [--linear-cfa] [--star-candidates] [--noise-model] [--frequency-cfa]\n"
      << "  " << executable
      << " develop <input-raw-or-dng> [--bias master.dng|--bias-dir dir] [--dark master.dng|--dark-dir dir] [--flat master.dng|--flat-dir dir] [--dark-excludes-bias] [--method bilinear-baseline|malvar-baseline|residual-interpolation|frequency-guided|inverse-refine] [-o output.tif|jpg] [--preview-stretch none|astro] [--jpeg-quality 1..100] [--inverse-iterations n] [--chroma-smoothness v] [--alias-suppression v] [--edge-sensitivity v] [--star-chroma-guard v] [--star-luma-threshold v] [--export-alias-risk map.tif|jpg] [--export-residual-map map.tif|jpg]\n"
      << "  " << executable
      << " calibrate <input-raw-or-dng> [--bias master.dng|--bias-dir dir] [--dark master.dng|--dark-dir dir] [--flat master.dng|--flat-dir dir] [--dark-excludes-bias] [-o output.tif|jpg] [--method inverse-refine]\n"
      << "  " << executable
      << " benchmark-debayer [--width n] [--height n] [--seed n] [--noise none|astro] [--export-prefix path]\n"
      << "  " << executable
      << " stack <input-raw-or-dng>... [--cfa-drizzle] [--auto-register] [--scale 2] [--offset dx,dy]...\n\n"
      << "Initial modes:\n"
      << "  faithful-astro   Conservative reconstruction from measured CFA samples\n"
      << "  star-preserve    Preserve stellar PSF shape and color discipline\n"
      << "  forensic         Emit residual, clipping, noise, and chroma-confidence maps\n";
}

inline void print_command_stub(std::ostream &out, std::string_view command) {
  out << "AstroCFA command '" << command
      << "' is scaffolded. The first implementation target is documented in "
         "docs/ROADMAP.md.\n";
}

inline bool is_known_command(std::string_view command) {
  return command == "inspect" || command == "develop" || command == "calibrate" ||
         command == "stack" || command == "benchmark-debayer";
}

} // namespace astrocfa
