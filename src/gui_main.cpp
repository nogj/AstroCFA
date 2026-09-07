#include "astrocfa/calibration.hpp"
#include "astrocfa/diagnostic_maps.hpp"
#include "astrocfa/demosaic.hpp"
#include "astrocfa/frequency_cfa.hpp"
#include "astrocfa/image_writer.hpp"
#include "astrocfa/output_transform.hpp"
#include "astrocfa/raw_color.hpp"
#include "astrocfa/raw_loader.hpp"
#include "astrocfa/noise_model.hpp"
#include "astrocfa/raw_inspector.hpp"
#include "astrocfa/star_detector.hpp"
#include "astrocfa/version.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPixmap>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QStatusBar>
#include <QStyle>
#include <QString>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

QPushButton *make_button(const QString &text, const QString &tooltip) {
  auto *button = new QPushButton(text);
  button->setToolTip(tooltip);
  button->setMinimumHeight(32);
  return button;
}

void append_log(QPlainTextEdit *log, const QString &message) {
  log->appendPlainText(message);
}

unsigned char to_display_u8(float value) {
  value = std::clamp(value, 0.0F, 1.0F);
  value = std::pow(value, 1.0F / 2.2F);
  return static_cast<unsigned char>(std::lround(value * 255.0F));
}

QImage to_qimage(const astrocfa::RgbImage &image) {
  QImage qimage(static_cast<int>(image.width()), static_cast<int>(image.height()),
                QImage::Format_RGB888);
  for(std::size_t y = 0; y < image.height(); ++y) {
    auto *row = qimage.scanLine(static_cast<int>(y));
    for(std::size_t x = 0; x < image.width(); ++x) {
      const astrocfa::RgbPixel pixel = image.pixel(x, y);
      row[x * 3U + 0U] = to_display_u8(pixel.r);
      row[x * 3U + 1U] = to_display_u8(pixel.g);
      row[x * 3U + 2U] = to_display_u8(pixel.b);
    }
  }
  return qimage;
}

std::string lowercase_extension(const std::string &path) {
  const std::size_t dot = path.find_last_of('.');
  if(dot == std::string::npos) {
    return {};
  }
  std::string extension = path.substr(dot + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
  return extension;
}

astrocfa::RgbImage gui_output_image(const astrocfa::RgbImage &linear,
                                    const std::string &path) {
  const std::string extension = lowercase_extension(path);
  if(extension == "jpg" || extension == "jpeg") {
    return astrocfa::make_astro_preview(linear);
  }
  return linear;
}

bool is_raw_like_path(const std::filesystem::path &path) {
  std::string extension = path.extension().string();
  std::transform(extension.begin(), extension.end(), extension.begin(),
                 [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
  return extension == ".dng" || extension == ".cr2" || extension == ".cr3" ||
         extension == ".nef" || extension == ".arw" || extension == ".raf" ||
         extension == ".orf" || extension == ".rw2" || extension == ".pef" ||
         extension == ".srw";
}

struct LoadedMaster {
  std::unique_ptr<astrocfa::LinearRawFrame> file;
  std::unique_ptr<astrocfa::MasterBuildResult> directory;

  [[nodiscard]] const astrocfa::CfaFrame *cfa() const {
    return directory ? &directory->cfa : file ? &file->cfa : nullptr;
  }
};

LoadedMaster load_master(const QString &selection) {
  LoadedMaster loaded;
  if(selection.isEmpty()) {
    return loaded;
  }

  const std::filesystem::path path(selection.toStdString());
  if(!std::filesystem::is_directory(path)) {
    loaded.file = std::make_unique<astrocfa::LinearRawFrame>(
        astrocfa::load_linear_cfa_file(path.string()));
    return loaded;
  }

  std::vector<std::filesystem::path> paths;
  for(const auto &entry : std::filesystem::directory_iterator(path)) {
    if(entry.is_regular_file() && is_raw_like_path(entry.path())) {
      paths.push_back(entry.path());
    }
  }
  std::sort(paths.begin(), paths.end());
  if(paths.empty()) {
    throw std::invalid_argument("No RAW/DNG files found in master directory: " +
                                path.string());
  }

  std::vector<astrocfa::CfaFrame> frames;
  frames.reserve(paths.size());
  for(const auto &frame_path : paths) {
    frames.push_back(astrocfa::load_linear_cfa_file(frame_path.string()).cfa);
  }
  loaded.directory = std::make_unique<astrocfa::MasterBuildResult>(
      astrocfa::build_master_cfa(frames));
  return loaded;
}

void write_master_report(const char *name, const LoadedMaster &master,
                         std::ostream &out) {
  if(master.directory) {
    const auto &stats = master.directory->stats;
    out << "  " << name << " master: " << stats.frames << " frames, mean "
        << std::fixed << std::setprecision(8) << stats.mean << ", "
        << stats.rejected_samples << " rejected samples\n";
  } else if(master.file) {
    out << "  " << name << " master: file\n";
  }
}

struct ReconstructionPreset {
  std::string name;
  astrocfa::DemosaicResult result;
};

ReconstructionPreset reconstruct_for_preset(const QString &preset,
                                             const astrocfa::CfaFrame &cfa) {
  const astrocfa::NoiseModel noise;
  if(preset == "faithful-astro") {
    return {"inverse-refine", astrocfa::reconstruct_inverse_refine(cfa, noise)};
  }
  if(preset == "star-preserve") {
    astrocfa::InverseRefinementOptions options;
    options.star_chroma_guard = 0.60;
    options.alias_suppression = 0.55;
    return {"inverse-refine (star-preserve)",
            astrocfa::reconstruct_inverse_refine(cfa, noise, options)};
  }
  if(preset == "forensic" || preset == "frequency-guided") {
    return {preset == "forensic" ? "frequency-guided (forensic)" : "frequency-guided",
            astrocfa::reconstruct_frequency_guided(cfa, noise)};
  }
  if(preset == "inverse-refine") {
    return {"inverse-refine", astrocfa::reconstruct_inverse_refine(cfa, noise)};
  }
  if(preset == "malvar-baseline") {
    return {"malvar-baseline", astrocfa::reconstruct_malvar_baseline(cfa, noise)};
  }
  if(preset == "residual-interpolation") {
    return {"residual-interpolation",
            astrocfa::reconstruct_residual_interpolation(cfa, noise)};
  }
  return {"bilinear-baseline", astrocfa::reconstruct_baseline(cfa, noise)};
}

} // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QApplication::setApplicationName("AstroCFA");
  QApplication::setApplicationVersion(astrocfa::version);

  QMainWindow window;
  window.setWindowTitle("AstroCFA");
  window.resize(1040, 900);

  auto *central = new QWidget;
  auto *root = new QVBoxLayout(central);
  root->setContentsMargins(16, 16, 16, 16);
  root->setSpacing(12);

  auto *title = new QLabel("AstroCFA");
  QFont title_font = title->font();
  title_font.setPointSize(22);
  title_font.setBold(true);
  title->setFont(title_font);

  auto *subtitle =
      new QLabel("CFA-aware, measurement-constrained RAW development for astrophotography");
  subtitle->setWordWrap(true);

  root->addWidget(title);
  root->addWidget(subtitle);

  auto *input_group = new QGroupBox("Input");
  auto *input_layout = new QGridLayout(input_group);
  auto *input_path = new QLineEdit;
  auto *output_path = new QLineEdit;
  input_path->setPlaceholderText("RAW, DNG, EXR, or calibrated frame");
  output_path->setPlaceholderText("Optional TIFF/JPEG output");
  auto *browse = make_button("Browse", "Choose the source frame or master");
  auto *browse_output = make_button("Output", "Choose TIFF or JPEG output");
  input_layout->addWidget(new QLabel("Source"), 0, 0);
  input_layout->addWidget(input_path, 0, 1);
  input_layout->addWidget(browse, 0, 2);
  input_layout->addWidget(new QLabel("Output"), 1, 0);
  input_layout->addWidget(output_path, 1, 1);
  input_layout->addWidget(browse_output, 1, 2);

  auto *calibration_group = new QGroupBox("Calibration masters");
  auto *calibration_layout = new QGridLayout(calibration_group);
  auto *bias_path = new QLineEdit;
  auto *dark_path = new QLineEdit;
  auto *flat_path = new QLineEdit;
  auto *dark_includes_bias = new QCheckBox("Dark includes bias");
  auto *cosmetic_correction = new QCheckBox("Cosmetic correction");
  dark_includes_bias->setChecked(true);
  cosmetic_correction->setChecked(true);
  dark_includes_bias->setToolTip(
      "Disable only when the master dark was built after bias subtraction");
  cosmetic_correction->setToolTip(
      "Detect defects from dark/flat masters and repair from same-phase CFA neighbors");
  const auto add_master_row = [&](int row, const QString &label, QLineEdit *path) {
    path->setPlaceholderText("Optional master file or RAW directory");
    path->setClearButtonEnabled(true);
    auto *choose_file = new QToolButton;
    choose_file->setIcon(window.style()->standardIcon(QStyle::SP_FileIcon));
    choose_file->setToolTip("Choose " + label.toLower() + " master file");
    auto *choose_dir = new QToolButton;
    choose_dir->setIcon(window.style()->standardIcon(QStyle::SP_DirOpenIcon));
    choose_dir->setToolTip("Build " + label.toLower() + " master from a directory");
    calibration_layout->addWidget(new QLabel(label), row, 0);
    calibration_layout->addWidget(path, row, 1);
    calibration_layout->addWidget(choose_file, row, 2);
    calibration_layout->addWidget(choose_dir, row, 3);
    QObject::connect(choose_file, &QToolButton::clicked, [&, path, label]() {
      const QString file = QFileDialog::getOpenFileName(
          &window, "Choose " + label.toLower() + " master", QString(),
          "RAW masters (*.dng *.cr2 *.cr3 *.nef *.arw *.raf *.orf *.rw2 *.pef *.srw);;All files (*)");
      if(!file.isEmpty()) {
        path->setText(file);
      }
    });
    QObject::connect(choose_dir, &QToolButton::clicked, [&, path, label]() {
      const QString directory = QFileDialog::getExistingDirectory(
          &window, "Choose " + label.toLower() + " frames directory");
      if(!directory.isEmpty()) {
        path->setText(directory);
      }
    });
  };
  add_master_row(0, "Bias", bias_path);
  add_master_row(1, "Dark", dark_path);
  add_master_row(2, "Flat", flat_path);
  calibration_layout->addWidget(dark_includes_bias, 3, 1);
  calibration_layout->addWidget(cosmetic_correction, 3, 2, 1, 2);

  auto *mode_group = new QGroupBox("Mode");
  auto *mode_layout = new QGridLayout(mode_group);
  auto *mode = new QComboBox;
  auto *linear_cfa = new QCheckBox("Linear CFA check");
  auto *star_candidates = new QCheckBox("Star candidates");
  auto *noise_model = new QCheckBox("Noise model");
  auto *frequency_cfa = new QCheckBox("Frequency CFA");
  auto *white_balance = new QComboBox;
  auto *output_space = new QComboBox;
  linear_cfa->setToolTip("Load normalized CFA samples without demosaicing");
  star_candidates->setToolTip("Detect bright candidates on a CFA-safe luminance proxy");
  noise_model->setToolTip("Show the initial Poisson-Gaussian noise model");
  frequency_cfa->setToolTip("Estimate Bayer carrier energy and alias risk before demosaicing");
  mode->addItem("faithful-astro");
  mode->addItem("star-preserve");
  mode->addItem("forensic");
  mode->addItem("inverse-refine");
  mode->addItem("frequency-guided");
  mode->addItem("malvar-baseline");
  mode->addItem("residual-interpolation");
  white_balance->addItem("auto");
  white_balance->addItem("as-shot");
  white_balance->addItem("daylight");
  white_balance->addItem("unity");
  white_balance->setToolTip("Select RAW metadata or unity sensor balance");
  output_space->addItem("auto");
  output_space->addItem("linear sRGB");
  output_space->addItem("camera RGB");
  output_space->setToolTip("Convert with the camera matrix or preserve sensor RGB");
  mode_layout->addWidget(new QLabel("Reconstruction"), 0, 0);
  mode_layout->addWidget(mode, 0, 1);
  mode_layout->addWidget(linear_cfa, 0, 2);
  mode_layout->addWidget(star_candidates, 0, 3);
  mode_layout->addWidget(noise_model, 0, 4);
  mode_layout->addWidget(frequency_cfa, 0, 5);
  mode_layout->addWidget(new QLabel("White balance"), 1, 0);
  mode_layout->addWidget(white_balance, 1, 1);
  mode_layout->addWidget(new QLabel("Output space"), 1, 2);
  mode_layout->addWidget(output_space, 1, 3);

  auto *actions = new QWidget;
  auto *actions_layout = new QHBoxLayout(actions);
  actions_layout->setContentsMargins(0, 0, 0, 0);
  auto *inspect = make_button("Inspect", "Read RAW/CFA metadata and diagnostics");
  auto *calibrate = make_button("Calibrate", "Apply bias, dark, flat, and cosmetic correction");
  auto *stack = make_button("Stack", "Register and integrate frames with optional CFA drizzle");
  auto *develop = make_button("Develop", "Reconstruct and export the selected image");
  actions_layout->addWidget(inspect);
  actions_layout->addWidget(calibrate);
  actions_layout->addWidget(stack);
  actions_layout->addWidget(develop);
  actions_layout->addStretch(1);

  auto *preview_group = new QGroupBox("Preview");
  auto *preview_layout = new QVBoxLayout(preview_group);
  auto *preview_toolbar = new QWidget;
  auto *preview_toolbar_layout = new QHBoxLayout(preview_toolbar);
  preview_toolbar_layout->setContentsMargins(0, 0, 0, 0);
  auto *overlay = new QComboBox;
  overlay->addItem("image");
  overlay->addItem("alias risk");
  overlay->addItem("residual");
  overlay->addItem("sensor defects");
  overlay->setToolTip(
      "Sensor defects: hot red, dead blue, mixed magenta, invalid master yellow");
  auto *zoom = new QComboBox;
  zoom->addItem("fit");
  zoom->addItem("100%");
  zoom->addItem("200%");
  auto *preview_status = new QLabel("No preview");
  preview_toolbar_layout->addWidget(new QLabel("Overlay"));
  preview_toolbar_layout->addWidget(overlay);
  preview_toolbar_layout->addWidget(new QLabel("Zoom"));
  preview_toolbar_layout->addWidget(zoom);
  preview_toolbar_layout->addStretch(1);
  preview_toolbar_layout->addWidget(preview_status);

  auto *image_label = new QLabel;
  image_label->setAlignment(Qt::AlignCenter);
  image_label->setMinimumSize(480, 240);
  image_label->setText("Develop an input to preview reconstruction and diagnostics");
  auto *scroll_area = new QScrollArea;
  scroll_area->setWidget(image_label);
  scroll_area->setWidgetResizable(true);
  preview_layout->addWidget(preview_toolbar);
  preview_layout->addWidget(scroll_area, 1);

  auto *log = new QPlainTextEdit;
  log->setReadOnly(true);
  log->setMinimumHeight(140);
  log->setPlaceholderText("Processing log");

  root->addWidget(input_group);
  root->addWidget(calibration_group);
  root->addWidget(mode_group);
  root->addWidget(actions);
  root->addWidget(preview_group, 2);
  root->addWidget(log, 1);

  window.setCentralWidget(central);
  window.statusBar()->showMessage("Ready");

  QObject::connect(browse, &QPushButton::clicked, [&]() {
    const QString file = QFileDialog::getOpenFileName(
        &window, "Open AstroCFA input", QString(),
        "Astro inputs (*.dng *.cr2 *.cr3 *.nef *.arw *.raf *.fits *.fit *.exr *.tif *.tiff);;All files (*)");
    if(!file.isEmpty()) {
      input_path->setText(file);
      append_log(log, "Selected: " + file);
    }
  });

  QObject::connect(browse_output, &QPushButton::clicked, [&]() {
    const QString file = QFileDialog::getSaveFileName(
        &window, "Save AstroCFA output", QString(),
        "TIFF image (*.tif *.tiff);;JPEG preview (*.jpg *.jpeg);;All files (*)");
    if(!file.isEmpty()) {
      output_path->setText(file);
      append_log(log, "Output: " + file);
    }
  });

  const auto require_input = [&]() {
    if(input_path->text().isEmpty()) {
      append_log(log, "Choose an input first.");
      window.statusBar()->showMessage("Input required");
      return false;
    }
    return true;
  };

  auto preview_image = std::make_shared<QImage>();
  auto alias_image = std::make_shared<QImage>();
  auto residual_image = std::make_shared<QImage>();
  auto defect_image = std::make_shared<QImage>();

  const auto update_preview = [&]() {
    const QImage *selected = nullptr;
    if(overlay->currentText() == "alias risk") {
      selected = alias_image.get();
    } else if(overlay->currentText() == "residual") {
      selected = residual_image.get();
    } else if(overlay->currentText() == "sensor defects") {
      selected = defect_image.get();
    } else {
      selected = preview_image.get();
    }

    if(selected == nullptr || selected->isNull()) {
      image_label->setPixmap(QPixmap());
      image_label->setText("Develop an input to preview reconstruction and diagnostics");
      preview_status->setText("No preview");
      return;
    }

    QPixmap pixmap = QPixmap::fromImage(*selected);
    if(zoom->currentText() == "fit") {
      const QSize target = scroll_area->viewport()->size();
      pixmap = pixmap.scaled(target, Qt::KeepAspectRatio, Qt::SmoothTransformation);
    } else if(zoom->currentText() == "200%") {
      pixmap = pixmap.scaled(selected->width() * 2, selected->height() * 2,
                             Qt::KeepAspectRatio, Qt::FastTransformation);
    }
    image_label->setText(QString());
    image_label->setPixmap(pixmap);
    image_label->resize(pixmap.size());
    preview_status->setText(QString("%1 x %2").arg(selected->width()).arg(selected->height()));
  };

  QObject::connect(overlay, &QComboBox::currentTextChanged, [&]() { update_preview(); });
  QObject::connect(zoom, &QComboBox::currentTextChanged, [&]() { update_preview(); });

  QObject::connect(inspect, &QPushButton::clicked, [&]() {
    if(!require_input()) {
      return;
    }

    append_log(log, "inspect: reading RAW/CFA metadata...");
    try {
      const auto inspection = astrocfa::inspect_raw_file(input_path->text().toStdString());
      std::ostringstream report;
      astrocfa::write_inspection_report(inspection, report);
      if(linear_cfa->isChecked() || star_candidates->isChecked() ||
         noise_model->isChecked() || frequency_cfa->isChecked()) {
        const auto frame = astrocfa::load_linear_cfa_file(input_path->text().toStdString());
        if(linear_cfa->isChecked()) {
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
        report << "\nLinear CFA load check:\n"
               << "  dimensions: " << frame.cfa.width() << " x " << frame.cfa.height()
               << "\n"
               << "  valid samples: " << valid << "\n"
               << "  clipped samples: " << clipped << "\n"
               << "  mean unclipped normalized signal: " << std::fixed
               << std::setprecision(6) << (usable > 0.0 ? sum / usable : 0.0) << "\n";
        }
        if(star_candidates->isChecked()) {
          const astrocfa::StarDetectionStats stars = astrocfa::detect_star_candidates(frame.cfa);
          report << "\nCFA-safe star candidate check:\n"
                 << "  background mean: " << std::fixed << std::setprecision(6)
                 << stars.background_mean << "\n"
                 << "  background sigma: " << stars.background_sigma << "\n"
                 << "  detection threshold: " << stars.threshold << "\n"
                 << "  candidates: " << stars.candidates << "\n"
                 << "  largest candidate area: " << stars.largest_area << " proxy samples\n"
                 << "  brightest proxy signal: " << stars.brightest << "\n";
        }
        if(noise_model->isChecked()) {
          const astrocfa::NoiseModel model;
          const auto sky = astrocfa::estimate_noise(0.01, model);
          const auto mid = astrocfa::estimate_noise(0.25, model);
          const auto bright = astrocfa::estimate_noise(0.75, model);
          report << "\nInitial Poisson-Gaussian noise model:\n"
                 << "  read noise: " << model.read_noise << " normalized units\n"
                 << "  shot noise scale: " << model.shot_noise_scale << "\n"
                 << "  sigma @ 1% signal: " << sky.sigma << "\n"
                 << "  sigma @ 25% signal: " << mid.sigma << "\n"
                 << "  sigma @ 75% signal: " << bright.sigma << "\n";
        }
        if(frequency_cfa->isChecked()) {
          const astrocfa::FrequencyCfaDiagnostics frequency =
              astrocfa::analyze_frequency_cfa(frame.cfa);
          const double high_risk_percent =
              frequency.total_tiles > 0
                  ? 100.0 * static_cast<double>(frequency.high_risk_tiles) /
                        static_cast<double>(frequency.total_tiles)
                  : 0.0;
          report << "\nCFA frequency diagnostics:\n"
                 << "  tile size: " << frequency.tile_size << "\n"
                 << "  tiles: " << frequency.total_tiles << "\n"
                 << "  high-risk tiles: " << frequency.high_risk_tiles << " ("
                 << std::fixed << std::setprecision(2) << high_risk_percent << "%)\n"
                 << "  mean AC energy: " << std::setprecision(8)
                 << frequency.mean_ac_energy << "\n"
                 << "  mean CFA carrier energy: " << frequency.mean_carrier_energy << "\n"
                 << "  mean alias risk: " << frequency.mean_alias_risk << "\n"
                 << "  max alias risk: " << frequency.max_alias_risk << "\n";
        }
      }
      append_log(log, QString::fromStdString(report.str()));
      window.statusBar()->showMessage("Inspection complete");
    } catch(const std::exception &error) {
      append_log(log, "inspect failed: " + QString::fromUtf8(error.what()));
      window.statusBar()->showMessage("Inspection failed");
    }
  });

  const auto run_stub = [&](const QString &command) {
    if(!require_input()) {
      return;
    }
    append_log(log, command + ": scaffolded for " + input_path->text() +
                        " using mode " + mode->currentText() + ".");
    window.statusBar()->showMessage(command + " scaffolded");
  };

  QObject::connect(stack, &QPushButton::clicked, [&]() { run_stub("stack"); });

  const auto run_reconstruction = [&](bool export_result) {
    if(!require_input()) {
      return;
    }

    const QString action = export_result ? "develop" : "calibrate";
    append_log(log, action + ": loading linear CFA and calibration masters...");
    QApplication::setOverrideCursor(Qt::WaitCursor);
    try {
      const auto frame = astrocfa::load_linear_cfa_file(input_path->text().toStdString());
      const LoadedMaster bias = load_master(bias_path->text());
      const LoadedMaster dark = load_master(dark_path->text());
      const LoadedMaster flat = load_master(flat_path->text());
      const astrocfa::CalibrationResult calibrated = astrocfa::calibrate_cfa(
          frame.cfa,
          astrocfa::CalibrationInputs{
              .bias = bias.cfa(),
              .dark = dark.cfa(),
              .flat = flat.cfa(),
              .options = astrocfa::CalibrationOptions{
                  .dark_includes_bias = dark_includes_bias->isChecked(),
                  .cosmetic = astrocfa::CosmeticCorrectionOptions{
                      .enabled = cosmetic_correction->isChecked(),
                  },
              },
          });
      const ReconstructionPreset reconstruction =
          reconstruct_for_preset(mode->currentText(), calibrated.cfa);
      const astrocfa::DemosaicResult &result = reconstruction.result;
      astrocfa::RawColorOptions color_options =
          astrocfa::automatic_raw_color_options(frame.inspection.color);
      if(white_balance->currentText() == "as-shot") {
        color_options.white_balance = astrocfa::WhiteBalanceMode::as_shot;
      } else if(white_balance->currentText() == "daylight") {
        color_options.white_balance = astrocfa::WhiteBalanceMode::daylight;
      } else if(white_balance->currentText() == "unity") {
        color_options.white_balance = astrocfa::WhiteBalanceMode::unity;
      }
      if(output_space->currentText() != "auto") {
        color_options.convert_to_srgb =
            output_space->currentText() == "linear sRGB";
      }
      const astrocfa::RawColorResult developed = astrocfa::apply_raw_color(
          result.image, frame.inspection.color, color_options);
      std::ostringstream report;
      const astrocfa::DemosaicQuality quality =
          astrocfa::analyze_demosaic_quality(result.image, calibrated.cfa);
      *preview_image = to_qimage(astrocfa::make_astro_preview(developed.image));
      *alias_image = to_qimage(astrocfa::make_frequency_alias_risk_map(calibrated.cfa));
      *residual_image = to_qimage(
          astrocfa::make_remosaic_residual_map(calibrated.cfa, result.image));
      *defect_image = to_qimage(astrocfa::make_sensor_defect_map(calibrated.defects));
      update_preview();
      report << "AstroCFA calibrated reconstruction\n"
             << "  input: " << input_path->text().toStdString() << "\n"
             << "  preset: " << mode->currentText().toStdString() << "\n"
             << "  method: " << reconstruction.name << "\n"
             << "  dimensions: " << calibrated.cfa.width() << " x "
             << calibrated.cfa.height() << "\n";
      write_master_report("bias", bias, report);
      write_master_report("dark", dark, report);
      write_master_report("flat", flat, report);
      report << "  calibration samples: " << calibrated.stats.samples << "\n"
             << "  invalid samples: " << calibrated.stats.invalid_samples << "\n"
             << "  clipped samples: " << calibrated.stats.clipped_samples << "\n"
             << "  flat floor samples: " << calibrated.stats.flat_floor_samples << "\n"
             << "  hot/dead pixels: " << calibrated.stats.hot_pixels << " / "
             << calibrated.stats.dead_pixels << "\n"
             << "  invalid master pixels: "
             << calibrated.stats.invalid_master_pixels << "\n"
             << "  repaired/unrepaired defects: " << calibrated.stats.repaired_pixels
             << " / " << calibrated.stats.unrepaired_pixels << "\n"
             << "  mean before/after: " << std::fixed << std::setprecision(8)
             << calibrated.stats.mean_before << " / " << calibrated.stats.mean_after << "\n"
             << "  mean bias subtracted: " << calibrated.stats.mean_bias_subtracted << "\n"
             << "  mean dark subtracted: " << calibrated.stats.mean_dark_subtracted << "\n"
             << "  remosaic residual samples: " << result.residual.samples << "\n"
             << "  remosaic residual MAE: "
             << result.residual.mean_absolute << "\n"
             << "  remosaic residual RMS: " << result.residual.root_mean_square << "\n"
             << "  remosaic residual max: " << result.residual.maximum_absolute << "\n"
             << "  reduced chi-square: "
             << result.noise_weighted_residual.reduced_chi_square << "\n"
             << "  mean chroma roughness: " << quality.mean_chroma_roughness
             << "\n"
             << "  mean interpolated chroma: " << quality.mean_interpolated_chroma
             << "\n"
             << "  white balance: " << developed.stats.white_balance[0] << ", "
             << developed.stats.white_balance[1] << ", "
             << developed.stats.white_balance[2] << "\n"
             << "  output color space: "
             << (developed.stats.converted_to_srgb ? "linear sRGB"
                                                   : "camera RGB")
             << "\n"
             << "  negative/out-of-range pixels: "
             << developed.stats.negative_pixels << " / "
             << developed.stats.over_range_pixels << "\n";
      if(export_result && !output_path->text().isEmpty()) {
        const std::string path = output_path->text().toStdString();
        astrocfa::write_rgb_image(gui_output_image(developed.image, path), path);
        report << "  output: " << output_path->text().toStdString() << "\n";
      } else if(export_result) {
        report << "  note: no output path selected; choose one to write TIFF/JPEG.\n";
      }
      append_log(log, QString::fromStdString(report.str()));
      window.statusBar()->showMessage(export_result ? "Develop complete"
                                                    : "Calibration preview complete");
    } catch(const std::exception &error) {
      append_log(log, action + " failed: " + QString::fromUtf8(error.what()));
      window.statusBar()->showMessage(action + " failed");
    }
    QApplication::restoreOverrideCursor();
  };

  QObject::connect(calibrate, &QPushButton::clicked,
                   [&]() { run_reconstruction(false); });
  QObject::connect(develop, &QPushButton::clicked,
                   [&]() { run_reconstruction(true); });

  window.show();
  return QApplication::exec();
}
