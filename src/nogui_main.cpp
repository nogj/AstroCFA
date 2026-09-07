#include "astrocfa/cli.hpp"
#include "astrocfa/calibration.hpp"
#include "astrocfa/diagnostic_maps.hpp"
#include "astrocfa/demosaic.hpp"
#include "astrocfa/drizzle.hpp"
#include "astrocfa/frequency_cfa.hpp"
#include "astrocfa/image_writer.hpp"
#include "astrocfa/noise_model.hpp"
#include "astrocfa/output_transform.hpp"
#include "astrocfa/raw_loader.hpp"
#include "astrocfa/raw_inspector.hpp"
#include "astrocfa/registration.hpp"
#include "astrocfa/star_detector.hpp"
#include "astrocfa/version.hpp"

#include <iomanip>
#include <cctype>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

astrocfa::SubpixelOffset parse_offset(const std::string &text) {
  const std::size_t comma = text.find(',');
  if(comma == std::string::npos) {
    throw std::invalid_argument("Offset must use dx,dy format");
  }

  return astrocfa::SubpixelOffset{
      .dx = std::stod(text.substr(0, comma)),
      .dy = std::stod(text.substr(comma + 1)),
  };
}

bool is_supported_demosaic_method(const std::string &method) {
  return method == "bilinear-baseline" || method == "malvar-baseline" ||
         method == "residual-interpolation" || method == "frequency-guided" ||
         method == "inverse-refine";
}

astrocfa::DemosaicResult reconstruct_with_method(const astrocfa::CfaFrame &cfa,
                                                 const std::string &method,
                                                 astrocfa::InverseRefinementOptions
                                                     inverse_options) {
  if(method == "bilinear-baseline") {
    return astrocfa::reconstruct_baseline(cfa, astrocfa::NoiseModel{});
  }
  if(method == "malvar-baseline") {
    return astrocfa::reconstruct_malvar_baseline(cfa, astrocfa::NoiseModel{});
  }
  if(method == "residual-interpolation") {
    return astrocfa::reconstruct_residual_interpolation(cfa, astrocfa::NoiseModel{});
  }
  if(method == "frequency-guided") {
    return astrocfa::reconstruct_frequency_guided(cfa, astrocfa::NoiseModel{});
  }
  if(method == "inverse-refine") {
    return astrocfa::reconstruct_inverse_refine(cfa, astrocfa::NoiseModel{},
                                                inverse_options);
  }
  throw std::invalid_argument("Unsupported demosaic method: " + method);
}

astrocfa::RgbImage output_image_for_path(const astrocfa::RgbImage &linear,
                                         const std::string &path,
                                         const std::string &preview_stretch) {
  const std::size_t dot = path.find_last_of('.');
  std::string extension = dot == std::string::npos ? std::string{} : path.substr(dot + 1);
  for(char &ch : extension) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }

  if((extension == "jpg" || extension == "jpeg") && preview_stretch == "astro") {
    return astrocfa::make_astro_preview(linear);
  }
  return linear;
}

struct CalibrationCliOptions {
  std::string bias_path;
  std::string dark_path;
  std::string flat_path;
  std::string bias_dir;
  std::string dark_dir;
  std::string flat_dir;
  astrocfa::CalibrationOptions options;
};

struct LoadedCalibration {
  std::unique_ptr<astrocfa::LinearRawFrame> bias;
  std::unique_ptr<astrocfa::LinearRawFrame> dark;
  std::unique_ptr<astrocfa::LinearRawFrame> flat;
  std::unique_ptr<astrocfa::MasterBuildResult> bias_master;
  std::unique_ptr<astrocfa::MasterBuildResult> dark_master;
  std::unique_ptr<astrocfa::MasterBuildResult> flat_master;
};

std::string lowercase(std::string value) {
  for(char &ch : value) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }
  return value;
}

bool is_raw_like_path(const std::filesystem::path &path) {
  const std::string extension = lowercase(path.extension().string());
  return extension == ".dng" || extension == ".cr2" || extension == ".cr3" ||
         extension == ".nef" || extension == ".arw" || extension == ".raf" ||
         extension == ".orf" || extension == ".rw2" || extension == ".pef" ||
         extension == ".srw";
}

std::vector<std::string> list_raw_like_files(const std::string &directory) {
  std::vector<std::string> paths;
  for(const std::filesystem::directory_entry &entry :
      std::filesystem::directory_iterator(directory)) {
    if(entry.is_regular_file() && is_raw_like_path(entry.path())) {
      paths.push_back(entry.path().string());
    }
  }
  std::sort(paths.begin(), paths.end());
  if(paths.empty()) {
    throw std::invalid_argument("No RAW/DNG files found in directory: " + directory);
  }
  return paths;
}

astrocfa::MasterBuildResult build_master_from_paths(const std::vector<std::string> &paths) {
  std::vector<astrocfa::CfaFrame> frames;
  frames.reserve(paths.size());
  for(const std::string &path : paths) {
    frames.push_back(astrocfa::load_linear_cfa_file(path).cfa);
  }
  return astrocfa::build_master_cfa(frames);
}

LoadedCalibration load_calibration_masters(const CalibrationCliOptions &options) {
  LoadedCalibration loaded;
  if(!options.bias_path.empty()) {
    loaded.bias = std::make_unique<astrocfa::LinearRawFrame>(
        astrocfa::load_linear_cfa_file(options.bias_path));
  }
  if(!options.dark_path.empty()) {
    loaded.dark = std::make_unique<astrocfa::LinearRawFrame>(
        astrocfa::load_linear_cfa_file(options.dark_path));
  }
  if(!options.flat_path.empty()) {
    loaded.flat = std::make_unique<astrocfa::LinearRawFrame>(
        astrocfa::load_linear_cfa_file(options.flat_path));
  }
  if(!options.bias_dir.empty()) {
    loaded.bias_master = std::make_unique<astrocfa::MasterBuildResult>(
        build_master_from_paths(list_raw_like_files(options.bias_dir)));
  }
  if(!options.dark_dir.empty()) {
    loaded.dark_master = std::make_unique<astrocfa::MasterBuildResult>(
        build_master_from_paths(list_raw_like_files(options.dark_dir)));
  }
  if(!options.flat_dir.empty()) {
    loaded.flat_master = std::make_unique<astrocfa::MasterBuildResult>(
        build_master_from_paths(list_raw_like_files(options.flat_dir)));
  }
  return loaded;
}

astrocfa::CalibrationResult apply_cli_calibration(
    const astrocfa::CfaFrame &light, const CalibrationCliOptions &options,
    const LoadedCalibration &loaded) {
  return astrocfa::calibrate_cfa(
      light, astrocfa::CalibrationInputs{
                 .bias = loaded.bias_master
                             ? &loaded.bias_master->cfa
                             : loaded.bias ? &loaded.bias->cfa : nullptr,
                 .dark = loaded.dark_master
                             ? &loaded.dark_master->cfa
                             : loaded.dark ? &loaded.dark->cfa : nullptr,
                 .flat = loaded.flat_master
                             ? &loaded.flat_master->cfa
                             : loaded.flat ? &loaded.flat->cfa : nullptr,
                 .options = options.options,
             });
}

void write_calibration_report(const astrocfa::CalibrationStats &stats, std::ostream &out) {
  out << "  calibration samples: " << stats.samples << "\n"
      << "  calibration invalid samples: " << stats.invalid_samples << "\n"
      << "  calibration clipped samples: " << stats.clipped_samples << "\n"
      << "  calibration flat floor samples: " << stats.flat_floor_samples << "\n"
      << "  calibration mean before: " << std::fixed << std::setprecision(8)
      << stats.mean_before << "\n"
      << "  calibration mean after: " << stats.mean_after << "\n"
      << "  mean bias subtracted: " << stats.mean_bias_subtracted << "\n"
      << "  mean dark subtracted: " << stats.mean_dark_subtracted << "\n"
      << "  flat phase means: " << stats.flat_phase_mean[0] << ", "
      << stats.flat_phase_mean[1] << ", " << stats.flat_phase_mean[2] << ", "
      << stats.flat_phase_mean[3] << "\n";
}

void write_master_report(const char *name, const astrocfa::MasterBuildResult *master,
                         std::ostream &out) {
  if(master == nullptr) {
    return;
  }
  out << "  " << name << " master frames: " << master->stats.frames << "\n"
      << "  " << name << " master mean: " << std::fixed << std::setprecision(8)
      << master->stats.mean << "\n"
      << "  " << name << " master invalid samples: "
      << master->stats.invalid_output_samples << "\n"
      << "  " << name << " master rejected samples: "
      << master->stats.rejected_samples << "\n";
}

void write_loaded_master_report(const LoadedCalibration &loaded, std::ostream &out) {
  write_master_report("bias", loaded.bias_master.get(), out);
  write_master_report("dark", loaded.dark_master.get(), out);
  write_master_report("flat", loaded.flat_master.get(), out);
}

void validate_calibration_cli_options(const CalibrationCliOptions &options) {
  if(!options.bias_path.empty() && !options.bias_dir.empty()) {
    throw std::invalid_argument("Use either --bias or --bias-dir, not both");
  }
  if(!options.dark_path.empty() && !options.dark_dir.empty()) {
    throw std::invalid_argument("Use either --dark or --dark-dir, not both");
  }
  if(!options.flat_path.empty() && !options.flat_dir.empty()) {
    throw std::invalid_argument("Use either --flat or --flat-dir, not both");
  }
}

} // namespace

int main(int argc, char **argv) {
  constexpr auto executable = "astrocfa-nogui";

  if(argc <= 1) {
    astrocfa::print_cli_help(std::cout, executable);
    return 0;
  }

  const std::string arg1 = argv[1];
  if(arg1 == "--help" || arg1 == "-h") {
    astrocfa::print_cli_help(std::cout, executable);
    return 0;
  }

  if(arg1 == "--version") {
    std::cout << astrocfa::version << "\n";
    return 0;
  }

  if(arg1 == "inspect") {
    if(argc < 3) {
      std::cerr << "inspect requires an input RAW/DNG path.\n";
      return 2;
    }

    try {
      const auto inspection = astrocfa::inspect_raw_file(argv[2]);
      astrocfa::write_inspection_report(inspection, std::cout);
      bool linear_cfa = false;
      bool star_candidates = false;
      bool noise_model = false;
      bool frequency_cfa = false;
      for(int i = 3; i < argc; ++i) {
        const std::string option = argv[i];
        linear_cfa = linear_cfa || option == "--linear-cfa";
        star_candidates = star_candidates || option == "--star-candidates";
        noise_model = noise_model || option == "--noise-model";
        frequency_cfa = frequency_cfa || option == "--frequency-cfa";
      }

      if(linear_cfa || star_candidates || noise_model || frequency_cfa) {
        const auto frame = astrocfa::load_linear_cfa_file(argv[2]);
        if(linear_cfa) {
          double sum = 0.0;
          std::size_t valid = 0;
          std::size_t clipped = 0;
          for(std::size_t y = 0; y < frame.cfa.height(); ++y) {
            for(std::size_t x = 0; x < frame.cfa.width(); ++x) {
              const auto sample = frame.cfa.sample_info(x, y);
              valid += sample.valid ? 1U : 0U;
              clipped += sample.clipped ? 1U : 0U;
              if(sample.valid && !sample.clipped) {
                sum += sample.value;
              }
            }
          }
          const double usable = static_cast<double>(valid - clipped);
          std::cout << "\nLinear CFA load check:\n"
                    << "  dimensions: " << frame.cfa.width() << " x " << frame.cfa.height()
                    << "\n"
                    << "  valid samples: " << valid << "\n"
                    << "  clipped samples: " << clipped << "\n"
                    << "  mean unclipped normalized signal: " << std::fixed
                    << std::setprecision(6) << (usable > 0.0 ? sum / usable : 0.0)
                    << "\n";
        }
        if(star_candidates) {
          const astrocfa::StarDetectionStats stars = astrocfa::detect_star_candidates(frame.cfa);
          std::cout << "\nCFA-safe star candidate check:\n"
                    << "  background mean: " << std::fixed << std::setprecision(6)
                    << stars.background_mean << "\n"
                    << "  background sigma: " << stars.background_sigma << "\n"
                    << "  detection threshold: " << stars.threshold << "\n"
                    << "  candidates: " << stars.candidates << "\n"
                    << "  largest candidate area: " << stars.largest_area
                    << " proxy samples\n"
                    << "  brightest proxy signal: " << stars.brightest << "\n";
        }
        if(noise_model) {
          const astrocfa::NoiseModel model;
          const auto sky = astrocfa::estimate_noise(0.01, model);
          const auto mid = astrocfa::estimate_noise(0.25, model);
          const auto bright = astrocfa::estimate_noise(0.75, model);
          std::cout << "\nInitial Poisson-Gaussian noise model:\n"
                    << "  read noise: " << model.read_noise << " normalized units\n"
                    << "  shot noise scale: " << model.shot_noise_scale << "\n"
                    << "  sigma @ 1% signal: " << sky.sigma << "\n"
                    << "  sigma @ 25% signal: " << mid.sigma << "\n"
                    << "  sigma @ 75% signal: " << bright.sigma << "\n";
        }
        if(frequency_cfa) {
          const astrocfa::FrequencyCfaDiagnostics frequency =
              astrocfa::analyze_frequency_cfa(frame.cfa);
          const double high_risk_percent =
              frequency.total_tiles > 0
                  ? 100.0 * static_cast<double>(frequency.high_risk_tiles) /
                        static_cast<double>(frequency.total_tiles)
                  : 0.0;
          std::cout << "\nCFA frequency diagnostics:\n"
                    << "  tile size: " << frequency.tile_size << "\n"
                    << "  tiles: " << frequency.total_tiles << "\n"
                    << "  high-risk tiles: " << frequency.high_risk_tiles << " ("
                    << std::fixed << std::setprecision(2) << high_risk_percent << "%)\n"
                    << "  mean AC energy: " << std::setprecision(8)
                    << frequency.mean_ac_energy << "\n"
                    << "  mean CFA carrier energy: " << frequency.mean_carrier_energy
                    << "\n"
                    << "  mean alias risk: " << frequency.mean_alias_risk << "\n"
                    << "  max alias risk: " << frequency.max_alias_risk << "\n";
        }
      }
    } catch(const std::exception &error) {
      std::cerr << "inspect failed: " << error.what() << "\n";
      return 1;
    }
    return 0;
  }

  if(arg1 == "develop") {
    if(argc < 3) {
      std::cerr << "develop requires an input RAW/DNG path.\n";
      return 2;
    }

    try {
      std::string method = "bilinear-baseline";
      std::string output_path;
      std::string alias_risk_path;
      std::string residual_map_path;
      std::string preview_stretch = "none";
      astrocfa::ImageWriteOptions write_options;
      astrocfa::InverseRefinementOptions inverse_options;
      CalibrationCliOptions calibration_options;
      for(int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if(arg == "--method" && i + 1 < argc) {
          method = argv[++i];
        } else if((arg == "-o" || arg == "--output") && i + 1 < argc) {
          output_path = argv[++i];
        } else if(arg == "--preview-stretch" && i + 1 < argc) {
          preview_stretch = argv[++i];
        } else if(arg == "--jpeg-quality" && i + 1 < argc) {
          write_options.jpeg_quality = std::stoi(argv[++i]);
        } else if(arg == "--inverse-iterations" && i + 1 < argc) {
          inverse_options.iterations = std::stoi(argv[++i]);
        } else if(arg == "--chroma-smoothness" && i + 1 < argc) {
          inverse_options.chroma_smoothness = std::stod(argv[++i]);
        } else if(arg == "--alias-suppression" && i + 1 < argc) {
          inverse_options.alias_suppression = std::stod(argv[++i]);
        } else if(arg == "--edge-sensitivity" && i + 1 < argc) {
          inverse_options.edge_sensitivity = std::stod(argv[++i]);
        } else if(arg == "--bias" && i + 1 < argc) {
          calibration_options.bias_path = argv[++i];
        } else if(arg == "--dark" && i + 1 < argc) {
          calibration_options.dark_path = argv[++i];
        } else if(arg == "--flat" && i + 1 < argc) {
          calibration_options.flat_path = argv[++i];
        } else if(arg == "--bias-dir" && i + 1 < argc) {
          calibration_options.bias_dir = argv[++i];
        } else if(arg == "--dark-dir" && i + 1 < argc) {
          calibration_options.dark_dir = argv[++i];
        } else if(arg == "--flat-dir" && i + 1 < argc) {
          calibration_options.flat_dir = argv[++i];
        } else if(arg == "--dark-excludes-bias") {
          calibration_options.options.dark_includes_bias = false;
        } else if(arg == "--export-alias-risk" && i + 1 < argc) {
          alias_risk_path = argv[++i];
        } else if(arg == "--export-residual-map" && i + 1 < argc) {
          residual_map_path = argv[++i];
        } else {
          throw std::invalid_argument("Unknown develop option: " + arg);
        }
      }

      if(!is_supported_demosaic_method(method)) {
        throw std::invalid_argument("Unsupported demosaic method: " + method);
      }
      if(preview_stretch != "none" && preview_stretch != "astro") {
        throw std::invalid_argument("Unsupported preview stretch: " + preview_stretch);
      }
      validate_calibration_cli_options(calibration_options);

      const auto frame = astrocfa::load_linear_cfa_file(argv[2]);
      const LoadedCalibration masters = load_calibration_masters(calibration_options);
      const astrocfa::CalibrationResult calibrated =
          apply_cli_calibration(frame.cfa, calibration_options, masters);
      astrocfa::DemosaicResult result =
          reconstruct_with_method(calibrated.cfa, method, inverse_options);
      std::cout << "AstroCFA reconstruction fidelity check\n"
                << "  input: " << argv[2] << "\n"
                << "  method: " << method << "\n"
                << "  dimensions: " << calibrated.cfa.width() << " x "
                << calibrated.cfa.height()
                << "\n"
                << "  inverse iterations: "
                << (method == "inverse-refine" ? inverse_options.iterations : 0)
                << "\n"
                << "  remosaic residual samples: " << result.residual.samples << "\n"
                << "  remosaic residual MAE: " << std::fixed << std::setprecision(8)
                << result.residual.mean_absolute << "\n"
                << "  remosaic residual RMS: " << result.residual.root_mean_square << "\n"
                << "  remosaic residual max: " << result.residual.maximum_absolute << "\n"
                << "  reduced chi-square: "
                << result.noise_weighted_residual.reduced_chi_square << "\n"
                << "  mean chroma roughness: ";
      const astrocfa::DemosaicQuality quality =
          astrocfa::analyze_demosaic_quality(result.image, calibrated.cfa);
      std::cout << quality.mean_chroma_roughness
                << "\n"
                << "  mean interpolated chroma: " << quality.mean_interpolated_chroma
                << "\n";
      if(!calibration_options.bias_path.empty() || !calibration_options.dark_path.empty() ||
         !calibration_options.flat_path.empty() || !calibration_options.bias_dir.empty() ||
         !calibration_options.dark_dir.empty() || !calibration_options.flat_dir.empty()) {
        write_loaded_master_report(masters, std::cout);
        write_calibration_report(calibrated.stats, std::cout);
      }
      if(!output_path.empty()) {
        astrocfa::write_rgb_image(
            output_image_for_path(result.image, output_path, preview_stretch),
            output_path, write_options);
        std::cout << "  output: " << output_path << "\n";
        if(preview_stretch == "astro") {
          std::cout << "  preview stretch: astro arcsinh for JPEG outputs\n";
        }
      } else {
        std::cout << "  note: no output path provided; use -o result.tif or -o preview.jpg.\n";
      }
      if(!alias_risk_path.empty()) {
        astrocfa::write_rgb_image(
            astrocfa::make_frequency_alias_risk_map(calibrated.cfa), alias_risk_path,
            write_options);
        std::cout << "  alias risk map: " << alias_risk_path << "\n";
      }
      if(!residual_map_path.empty()) {
        astrocfa::write_rgb_image(
            astrocfa::make_remosaic_residual_map(calibrated.cfa, result.image),
            residual_map_path, write_options);
        std::cout << "  remosaic residual map: " << residual_map_path << "\n";
      }
    } catch(const std::exception &error) {
      std::cerr << "develop failed: " << error.what() << "\n";
      return 1;
    }
    return 0;
  }

  if(arg1 == "calibrate") {
    if(argc < 3) {
      std::cerr << "calibrate requires an input RAW/DNG path.\n";
      return 2;
    }

    try {
      std::string method = "inverse-refine";
      std::string output_path;
      std::string preview_stretch = "none";
      astrocfa::ImageWriteOptions write_options;
      astrocfa::InverseRefinementOptions inverse_options;
      CalibrationCliOptions calibration_options;

      for(int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if(arg == "--method" && i + 1 < argc) {
          method = argv[++i];
        } else if((arg == "-o" || arg == "--output") && i + 1 < argc) {
          output_path = argv[++i];
        } else if(arg == "--preview-stretch" && i + 1 < argc) {
          preview_stretch = argv[++i];
        } else if(arg == "--jpeg-quality" && i + 1 < argc) {
          write_options.jpeg_quality = std::stoi(argv[++i]);
        } else if(arg == "--inverse-iterations" && i + 1 < argc) {
          inverse_options.iterations = std::stoi(argv[++i]);
        } else if(arg == "--chroma-smoothness" && i + 1 < argc) {
          inverse_options.chroma_smoothness = std::stod(argv[++i]);
        } else if(arg == "--alias-suppression" && i + 1 < argc) {
          inverse_options.alias_suppression = std::stod(argv[++i]);
        } else if(arg == "--edge-sensitivity" && i + 1 < argc) {
          inverse_options.edge_sensitivity = std::stod(argv[++i]);
        } else if(arg == "--bias" && i + 1 < argc) {
          calibration_options.bias_path = argv[++i];
        } else if(arg == "--dark" && i + 1 < argc) {
          calibration_options.dark_path = argv[++i];
        } else if(arg == "--flat" && i + 1 < argc) {
          calibration_options.flat_path = argv[++i];
        } else if(arg == "--bias-dir" && i + 1 < argc) {
          calibration_options.bias_dir = argv[++i];
        } else if(arg == "--dark-dir" && i + 1 < argc) {
          calibration_options.dark_dir = argv[++i];
        } else if(arg == "--flat-dir" && i + 1 < argc) {
          calibration_options.flat_dir = argv[++i];
        } else if(arg == "--dark-excludes-bias") {
          calibration_options.options.dark_includes_bias = false;
        } else {
          throw std::invalid_argument("Unknown calibrate option: " + arg);
        }
      }

      if(!is_supported_demosaic_method(method)) {
        throw std::invalid_argument("Unsupported demosaic method: " + method);
      }
      if(preview_stretch != "none" && preview_stretch != "astro") {
        throw std::invalid_argument("Unsupported preview stretch: " + preview_stretch);
      }
      validate_calibration_cli_options(calibration_options);

      const auto frame = astrocfa::load_linear_cfa_file(argv[2]);
      const LoadedCalibration masters = load_calibration_masters(calibration_options);
      const astrocfa::CalibrationResult calibrated =
          apply_cli_calibration(frame.cfa, calibration_options, masters);
      const astrocfa::DemosaicResult result =
          reconstruct_with_method(calibrated.cfa, method, inverse_options);
      const astrocfa::DemosaicQuality quality =
          astrocfa::analyze_demosaic_quality(result.image, calibrated.cfa);

      std::cout << "AstroCFA calibration and reconstruction\n"
                << "  input: " << argv[2] << "\n"
                << "  method: " << method << "\n";
      write_loaded_master_report(masters, std::cout);
      write_calibration_report(calibrated.stats, std::cout);
      std::cout << "  remosaic residual MAE: " << std::fixed << std::setprecision(8)
                << result.residual.mean_absolute << "\n"
                << "  remosaic residual RMS: " << result.residual.root_mean_square << "\n"
                << "  mean chroma roughness: " << quality.mean_chroma_roughness << "\n"
                << "  mean interpolated chroma: " << quality.mean_interpolated_chroma
                << "\n";
      if(!output_path.empty()) {
        astrocfa::write_rgb_image(
            output_image_for_path(result.image, output_path, preview_stretch),
            output_path, write_options);
        std::cout << "  output: " << output_path << "\n";
      } else {
        std::cout << "  note: no output path provided; use -o calibrated.tif or -o preview.jpg.\n";
      }
    } catch(const std::exception &error) {
      std::cerr << "calibrate failed: " << error.what() << "\n";
      return 1;
    }
    return 0;
  }

  if(arg1 == "stack") {
    std::vector<std::string> inputs;
    std::vector<astrocfa::SubpixelOffset> offsets;
    bool cfa_drizzle = false;
    bool auto_register = false;
    std::size_t scale = 2;

    for(int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if(arg == "--cfa-drizzle") {
        cfa_drizzle = true;
      } else if(arg == "--auto-register") {
        auto_register = true;
      } else if(arg == "--scale" && i + 1 < argc) {
        scale = static_cast<std::size_t>(std::stoul(argv[++i]));
      } else if(arg == "--offset" && i + 1 < argc) {
        offsets.push_back(parse_offset(argv[++i]));
      } else {
        inputs.push_back(arg);
      }
    }

    if(inputs.empty()) {
      std::cerr << "stack requires at least one input RAW/DNG path.\n";
      return 2;
    }
    if(!cfa_drizzle) {
      std::cout << "Only --cfa-drizzle stack mode is implemented at this stage.\n";
      return 0;
    }

    try {
      const auto first = astrocfa::load_linear_cfa_file(inputs.front());
      std::vector<astrocfa::SubpixelOffset> effective_offsets(inputs.size(),
                                                              astrocfa::SubpixelOffset{});
      std::vector<astrocfa::RegistrationResult> registrations(inputs.size());
      for(std::size_t i = 0; i < offsets.size() && i < effective_offsets.size(); ++i) {
        effective_offsets[i] = offsets[i];
      }

      astrocfa::CfaDrizzleAccumulator accumulator(
          first.cfa.width(), first.cfa.height(), first.cfa.pattern(),
          astrocfa::DrizzleOptions{.scale = scale, .drop_shrink = 0.7});
      accumulator.add_frame(first.cfa, effective_offsets[0]);

      for(std::size_t i = 1; i < inputs.size(); ++i) {
        const auto frame = astrocfa::load_linear_cfa_file(inputs[i]);
        if(auto_register && i >= offsets.size()) {
          registrations[i] =
              astrocfa::estimate_integer_star_translation(first.cfa, frame.cfa);
          effective_offsets[i] = registrations[i].offset;
        }
        accumulator.add_frame(frame.cfa, effective_offsets[i]);
      }

      const auto print_coverage = [](const char *name,
                                     const astrocfa::DrizzleCoverageStats &coverage) {
        const double percent =
            coverage.total_samples > 0
                ? 100.0 * static_cast<double>(coverage.covered_samples) /
                      static_cast<double>(coverage.total_samples)
                : 0.0;
        std::cout << "  " << name << ": " << coverage.covered_samples << "/"
                  << coverage.total_samples << " (" << std::fixed << std::setprecision(2)
                  << percent << "%), mean weight " << coverage.mean_weight
                  << ", max weight " << coverage.max_weight << "\n";
      };

      std::cout << "AstroCFA CFA drizzle accumulation\n"
                << "  frames: " << inputs.size() << "\n"
                << "  input dimensions: " << first.cfa.width() << " x " << first.cfa.height()
                << "\n"
                << "  output dimensions: " << accumulator.output_width() << " x "
                << accumulator.output_height() << "\n"
                << "  scale: " << scale << "\n"
                << "  offsets: " << offsets.size() << " provided"
                << (auto_register ? ", missing values estimated by CFA-safe star registration\n"
                                  : ", missing values default to 0,0\n");
      if(auto_register) {
        for(std::size_t i = 1; i < registrations.size(); ++i) {
          std::cout << "  registration frame " << i << ": dx="
                    << effective_offsets[i].dx << " dy=" << effective_offsets[i].dy
                    << " score=" << registrations[i].score
                    << " matched=" << registrations[i].matched_samples << "\n";
        }
      }
      std::cout
                << "  coverage by phase:\n";
      print_coverage("R ", accumulator.coverage(astrocfa::CfaColor::red));
      print_coverage("G1", accumulator.coverage(astrocfa::CfaColor::green1));
      print_coverage("B ", accumulator.coverage(astrocfa::CfaColor::blue));
      print_coverage("G2", accumulator.coverage(astrocfa::CfaColor::green2));
      std::cout << "  note: no image export yet; this validates phase-separated "
                   "drizzle accumulation.\n";
    } catch(const std::exception &error) {
      std::cerr << "stack failed: " << error.what() << "\n";
      return 1;
    }
    return 0;
  }

  if(astrocfa::is_known_command(arg1)) {
    astrocfa::print_command_stub(std::cout, arg1);
    return 0;
  }

  std::cerr << "Unknown command: " << arg1 << "\n\n";
  astrocfa::print_cli_help(std::cerr, executable);
  return 2;
}
