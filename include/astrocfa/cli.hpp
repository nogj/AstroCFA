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
      << " inspect <input-raw-or-dng> [--linear-cfa] [--star-candidates] [--estimate-psf] [--noise-model] [--frequency-cfa]\n"
      << "  " << executable
      << " develop <input-raw-or-dng> [--bias master.dng|--bias-dir dir] [--dark master.dng|--dark-dir dir] [--flat master.dng|--flat-dir dir] [--dark-excludes-bias] [--no-cosmetic-correction] [--background off|gradient|neutral] [--background-tile pixels] [--background-degree 0..2] [--white-balance auto|as-shot|daylight|unity] [--wb-multipliers r,g,b] [--output-space auto|srgb|camera] [--tone none|linear|arcsinh|ghs] [--exposure-ev v] [--black-point v --white-point v] [--stretch-factor v] [--local-intensity v] [--symmetry-point v] [--protect-shadows v] [--protect-highlights v] [--saturation v] [-o output.tif|jpg] [--preview-stretch none|astro] [--jpeg-quality 1..100] [--export-background map.tif|jpg] [--export-alias-risk map.tif|jpg] [--export-residual-map map.tif|jpg] [--export-defect-map map.tif|jpg]\n"
      << "  " << executable
      << " calibrate <input-raw-or-dng> [--bias master.dng|--bias-dir dir] [--dark master.dng|--dark-dir dir] [--flat master.dng|--flat-dir dir] [--dark-excludes-bias] [--no-cosmetic-correction] [--export-defect-map map.tif|jpg] [-o output.tif|jpg] [--preview-stretch none|astro] [--jpeg-quality 1..100]\n"
      << "  " << executable
      << " benchmark-debayer [--width n] [--height n] [--seed n] [--noise none|astro] [--stars n] [--common-star-sigma v] [--moffat-beta v] [--psf-ellipticity v] [--psf-angle radians] [--chromatic-psf-shift v] [--chromatic-psf-scale v] [--star-flux-scale v] [--star-phase diverse|redundant] [--export-prefix path] [--external name=linear-rgb.tif]... [--external-srgb name=srgb-rgb.tif]... [--csv results.csv]\n"
      << "  " << executable
      << " benchmark-joint [--width n] [--height n] [--frames n] [--iterations n] [--seed n] [--noise none|astro] [--seeing fixed|variable] [--transients n] [--luma-smoothness v] [--chroma-smoothness v] [--export-prefix path]\n"
      << "  " << executable
      << " stack <input-raw-or-dng>... [--joint-reconstruct] [--cfa-drizzle] [--auto-register] [--auto-psf|--psf-sigma pixels]... [--auto-psf-strength 0..1] [--scale 1|2] [--offset dx,dy]... [--iterations n] [--no-discrepancy-stop] [--minimum-iterations n] [--discrepancy-target v] [--luma-smoothness v] [--chroma-smoothness v] [--huber-sigma v] [--bias master.dng|--bias-dir dir] [--dark master.dng|--dark-dir dir] [--flat master.dng|--flat-dir dir] [-o output.tif|jpg] [--export-confidence map.tif|jpg]\n\n"
      << "Reconstruction:\n"
      << "  develop and calibrate automatically select the best validated AstroCFA\n"
      << "  point-source model (radial, achromatic ePSF, or chromatic ePSF).\n";
}

inline void print_command_stub(std::ostream &out, std::string_view command) {
  out << "AstroCFA command '" << command
      << "' is scaffolded. The first implementation target is documented in "
         "docs/ROADMAP.md.\n";
}

inline bool is_known_command(std::string_view command) {
  return command == "inspect" || command == "develop" || command == "calibrate" ||
         command == "stack" || command == "benchmark-debayer" ||
         command == "benchmark-joint";
}

} // namespace astrocfa
