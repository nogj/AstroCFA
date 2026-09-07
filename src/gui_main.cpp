#include "astrocfa/diagnostic_maps.hpp"
#include "astrocfa/demosaic.hpp"
#include "astrocfa/frequency_cfa.hpp"
#include "astrocfa/image_writer.hpp"
#include "astrocfa/output_transform.hpp"
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
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

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

} // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QApplication::setApplicationName("AstroCFA");
  QApplication::setApplicationVersion(astrocfa::version);

  QMainWindow window;
  window.setWindowTitle("AstroCFA");
  window.resize(920, 620);

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

  auto *mode_group = new QGroupBox("Mode");
  auto *mode_layout = new QGridLayout(mode_group);
  auto *mode = new QComboBox;
  auto *linear_cfa = new QCheckBox("Linear CFA check");
  auto *star_candidates = new QCheckBox("Star candidates");
  auto *noise_model = new QCheckBox("Noise model");
  auto *frequency_cfa = new QCheckBox("Frequency CFA");
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
  mode_layout->addWidget(new QLabel("Reconstruction"), 0, 0);
  mode_layout->addWidget(mode, 0, 1);
  mode_layout->addWidget(linear_cfa, 1, 1);
  mode_layout->addWidget(star_candidates, 2, 1);
  mode_layout->addWidget(noise_model, 3, 1);
  mode_layout->addWidget(frequency_cfa, 4, 1);

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
  image_label->setMinimumSize(480, 300);
  image_label->setText("Develop an input to preview reconstruction and diagnostics");
  auto *scroll_area = new QScrollArea;
  scroll_area->setWidget(image_label);
  scroll_area->setWidgetResizable(true);
  preview_layout->addWidget(preview_toolbar);
  preview_layout->addWidget(scroll_area, 1);

  auto *log = new QPlainTextEdit;
  log->setReadOnly(true);
  log->setMinimumHeight(220);
  log->setPlaceholderText("Processing log");

  root->addWidget(input_group);
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

  const auto update_preview = [&]() {
    const QImage *selected = nullptr;
    if(overlay->currentText() == "alias risk") {
      selected = alias_image.get();
    } else if(overlay->currentText() == "residual") {
      selected = residual_image.get();
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

  QObject::connect(calibrate, &QPushButton::clicked, [&]() { run_stub("calibrate"); });
  QObject::connect(stack, &QPushButton::clicked, [&]() { run_stub("stack"); });
  QObject::connect(develop, &QPushButton::clicked, [&]() {
    if(!require_input()) {
      return;
    }

    append_log(log, "develop: running baseline reconstruction fidelity check...");
    try {
      const auto frame = astrocfa::load_linear_cfa_file(input_path->text().toStdString());
      const bool use_inverse_refine = mode->currentText() == "inverse-refine";
      const bool use_frequency_guided = mode->currentText() == "frequency-guided";
      const bool use_malvar = mode->currentText() == "malvar-baseline";
      const bool use_ri = mode->currentText() == "residual-interpolation";
      const astrocfa::DemosaicResult result =
          use_inverse_refine
              ? astrocfa::reconstruct_inverse_refine(frame.cfa, astrocfa::NoiseModel{})
          : use_frequency_guided
              ? astrocfa::reconstruct_frequency_guided(frame.cfa, astrocfa::NoiseModel{})
              : use_malvar
                    ? astrocfa::reconstruct_malvar_baseline(frame.cfa, astrocfa::NoiseModel{})
              : use_ri ? astrocfa::reconstruct_residual_interpolation(frame.cfa,
                                                                      astrocfa::NoiseModel{})
                       : astrocfa::reconstruct_baseline(frame.cfa, astrocfa::NoiseModel{});
      std::ostringstream report;
      const astrocfa::DemosaicQuality quality =
          astrocfa::analyze_demosaic_quality(result.image, frame.cfa);
      *preview_image = to_qimage(astrocfa::make_astro_preview(result.image));
      *alias_image = to_qimage(astrocfa::make_frequency_alias_risk_map(frame.cfa));
      *residual_image = to_qimage(
          astrocfa::make_remosaic_residual_map(frame.cfa, result.image));
      update_preview();
      report << "AstroCFA reconstruction fidelity check\n"
             << "  input: " << input_path->text().toStdString() << "\n"
             << "  method: "
             << (use_inverse_refine ? "inverse-refine"
                 : use_frequency_guided ? "frequency-guided"
                                      : use_malvar   ? "malvar-baseline"
                                        : use_ri     ? "residual-interpolation"
                                                     : "bilinear-baseline")
             << "\n"
             << "  dimensions: " << frame.cfa.width() << " x " << frame.cfa.height() << "\n"
             << "  remosaic residual samples: " << result.residual.samples << "\n"
             << "  remosaic residual MAE: " << std::fixed << std::setprecision(8)
             << result.residual.mean_absolute << "\n"
             << "  remosaic residual RMS: " << result.residual.root_mean_square << "\n"
             << "  remosaic residual max: " << result.residual.maximum_absolute << "\n"
             << "  reduced chi-square: "
             << result.noise_weighted_residual.reduced_chi_square << "\n"
             << "  mean chroma roughness: " << quality.mean_chroma_roughness
             << "\n"
             << "  mean interpolated chroma: " << quality.mean_interpolated_chroma
             << "\n";
      if(!output_path->text().isEmpty()) {
        const std::string path = output_path->text().toStdString();
        astrocfa::write_rgb_image(gui_output_image(result.image, path), path);
        report << "  output: " << output_path->text().toStdString() << "\n";
      } else {
        report << "  note: no output path selected; choose one to write TIFF/JPEG.\n";
      }
      append_log(log, QString::fromStdString(report.str()));
      window.statusBar()->showMessage("Develop fidelity check complete");
    } catch(const std::exception &error) {
      append_log(log, "develop failed: " + QString::fromUtf8(error.what()));
      window.statusBar()->showMessage("Develop failed");
    }
  });

  window.show();
  return QApplication::exec();
}
