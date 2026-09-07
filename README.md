# AstroCFA

AstroCFA is a purpose-specific astrophotography engine for CFA-aware,
measurement-constrained RAW development.

It is not a general photo editor. Its goal is to extract the maximum defensible
astronomical signal from Bayer/X-Trans-like sensor data without neural networks,
while preserving traceability back to the measured samples.

```text
calibrated RAW lights
  -> CFA-aware registration and integration
  -> inverse demosaicing / reconstruction
  -> PSF-aware deconvolution
  -> photometric color calibration
  -> linear and stretched outputs
```

## Technical Thesis

AstroCFA treats demosaicing and development as constrained inverse problems.
The core question is:

```text
Which scene-linear RGB image could have produced these CFA measurements,
under the camera and noise model, while adding the least unsupported detail?
```

For deep-sky work this matters because tiny stars, weak nebulosity, dithering,
shot noise, hot pixels, clipping, and Bayer aliasing are not ordinary image
editing problems.

## Differentiators

- CFA drizzle and phase-aware integration before RGB reconstruction.
- Demosaicing refined by an explicit remosaicing residual.
- Poisson-Gaussian sensor noise modeling.
- Star-aware chroma discipline around undersampled and saturated stars.
- Explicit RAW white balance and camera-to-linear-sRGB color development.
- Classical PSF deconvolution, with no learned priors.
- Debug maps for clipping, residual, chroma confidence, aliasing, and noise.
- Reproducible command-line processing with sidecar metadata.

## Why AstroCFA Exists

Many astrophotography pipelines debayer early, then align and stack RGB frames.
AstroCFA's intended path is different:

```text
calibrate CFA lights
  -> register CFA-safe luminance proxies
  -> integrate measured CFA samples, optionally with drizzle
  -> reconstruct RGB only after the stack has accumulated real samples
```

For dithered deep-sky data this is the critical distinction. Dithered frames can
turn Bayer ambiguity into additional measured information. AstroCFA should use
that information before interpolation, then report where uncertainty remains.

## Initial App Shape

`astrocfa` is the Qt desktop application. It starts with a minimal interface for
selecting an input, loading bias/dark/flat masters, choosing a reconstruction
mode, previewing the reconstruction, switching diagnostic overlays, and exporting
TIFF/JPEG output. A calibration selector accepts either one master RAW or a
directory of RAW frames; directories are combined into a robust in-memory master.
The `faithful-astro` preset uses inverse refinement, `star-preserve` increases its
compact-star chroma guard, and `forensic` uses the frequency-guided path. White
balance and output-space selectors keep sensor RGB, as-shot/daylight metadata,
and linear sRGB conversion explicit.

`astrocfa-nogui` is the command-line executable for scripts, batch processing,
and reproducible long-running workflows.

## Initial CLI Shape

```bash
astrocfa-nogui inspect input.dng
astrocfa-nogui inspect input.dng --linear-cfa --star-candidates --estimate-psf \
  --noise-model --frequency-cfa
astrocfa-nogui benchmark-debayer
astrocfa-nogui benchmark-debayer --width 192 --height 128 --seed 7 --export-prefix bench/astro
astrocfa-nogui benchmark-joint --frames 4 --noise astro --seeing fixed
astrocfa-nogui benchmark-joint --seeing variable --transients 4 --export-prefix bench/joint
astrocfa-nogui calibrate light.dng --dark master-dark.dng --flat master-flat.dng \
  --method inverse-refine -o calibrated-preview.jpg --preview-stretch astro
astrocfa-nogui calibrate light.dng --dark-dir darks/ --flat-dir flats/ \
  --method inverse-refine -o calibrated-preview.jpg --preview-stretch astro
astrocfa-nogui stack light1.dng light2.dng light3.dng --cfa-drizzle --scale 2 \
  --offset 0,0 --offset 0.42,-0.18 --offset -0.31,0.27
astrocfa-nogui stack light1.dng light2.dng light3.dng --cfa-drizzle --auto-register
astrocfa-nogui stack light1.dng light2.dng light3.dng --joint-reconstruct \
  --offset 0,0 --offset 0.42,-0.18 --offset -0.31,0.27 --iterations 6 \
  -o joint-linear.tif --export-confidence joint-confidence.tif
astrocfa-nogui stack light1.dng light2.dng light3.dng --joint-reconstruct \
  --psf-sigma 0.72 --psf-sigma 0.91 --psf-sigma 0.68 -o joint-psf.tif
astrocfa-nogui stack light1.dng light2.dng light3.dng --joint-reconstruct \
  --auto-register --auto-psf -o joint-auto-psf.tif
astrocfa-nogui develop input.dng --method bilinear-baseline -o baseline.tif
astrocfa-nogui develop input.dng --method malvar-baseline -o malvar.tif
astrocfa-nogui develop input.dng --method residual-interpolation -o ri.tif
astrocfa-nogui develop input.dng --method frequency-guided -o preview.jpg
astrocfa-nogui develop input.dng --method inverse-refine -o preview.jpg \
  --preview-stretch astro --export-alias-risk alias.tif --export-residual-map residual.tif
astrocfa-nogui develop input.dng --method inverse-refine --white-balance daylight \
  --output-space srgb -o developed-linear-srgb.tif
astrocfa-nogui develop input.dng --method inverse-refine \
  --wb-multipliers 2.1,1.0,1.4 --output-space camera -o developed-camera-rgb.tif
astrocfa-nogui develop light.dng --bias master-bias.dng --dark master-dark.dng \
  --flat master-flat.dng --method inverse-refine -o calibrated.tif \
  --export-defect-map sensor-defects.tif
```

`develop` writes scene-linear 16-bit RGB TIFF output or an 8-bit JPEG preview.
Its default color mode uses as-shot WB and the camera-to-sRGB matrix when the RAW
provides them, with explicit daylight, unity, custom-multiplier, and camera-RGB
alternatives. JPEG output can use an astro-oriented arcsinh preview stretch while
TIFF remains linear. The `inverse-refine` method starts from frequency-guided residual
interpolation and iteratively regularizes chroma as an edge-aware field anchored
to measured CFA samples. It preserves the measured channel at every pixel, so
the remosaicing residual remains an explicit accountability check rather than a
decorative metric.

Diagnostic maps can be exported alongside the image. The first maps are Bayer
alias-risk and remosaicing residual. These are intentionally part of the core
workflow: AstroCFA should show which parts of the result are measured, inferred,
or risky.

`benchmark-debayer` generates a deterministic synthetic astro scene with known
RGB truth, mosaics it into CFA data, and compares every non-neural demosaic
candidate. It reports RGB error, chroma error, false star color, star luminance
error, aperture-flux error, FWHM error, elongation error, and remosaicing
residual. This is the early guardrail for keeping
AstroCFA's reconstruction work measurable rather than merely aesthetic.

`benchmark-joint` compares demosaic-each-light-then-average with phase-aware CFA
initialization, non-robust joint inversion, and robust joint inversion. It
synthesizes known dithers, Poisson-Gaussian noise, optional per-frame seeing
changes, and transient samples. The variable-seeing mode includes direct
ablations for ignoring the PSF, estimating it from CFA samples, and supplying
each light's known Gaussian PSF to the forward/adjoint operator.
The tuned priors remain exposed as `--luma-smoothness` and
`--chroma-smoothness` for reproducible ablations.

The inverse refinement path includes an optional star chroma guard for compact
PSF-like highlights. It is designed to reduce false magenta/green star cores
from unsupported Bayer chroma while preserving measured CFA samples exactly.

Calibration includes CFA-phase-aware cosmetic correction when a dark or flat
master is present. Hot pixels are local robust outliers in the dark; dead pixels
are local response deficits in the flat. Detected samples are replaced by the
median of valid, non-defective neighbors from the same Bayer phase before RGB
reconstruction. The map remains available as a Qt overlay and through
`--export-defect-map`; `--no-cosmetic-correction` preserves the uncorrected path.

`stack --joint-reconstruct` solves one RGB scene directly from all calibrated
CFA measurements and their supplied offsets. A phase-separated splat initializes
the scene, then robust noise-weighted backprojection reduces the residual against
every original light. Edge-aware chroma and green-luminance regularization keep
single-frame direct samples exact; with multiple noisy lights, measurements are
weighted fidelity constraints so the common estimate can denoise. The solver
models translation, bilinear sampling, and an optional Gaussian PSF per light.
Use repeated `--psf-sigma` values in sensor pixels (`FWHM / 2.35482`) when PSF
estimates are available. `--auto-psf` instead detects isolated unsaturated stars
on a CFA-safe luminance proxy, fits their profiles against the original mosaic
samples, and models each light's blur relative to the sharpest frame. Its
default strength is deliberately conservative and can be reproduced with
`--auto-psf-strength 0.75`.

Joint reconstruction stops by default when the robust reduced chi-square reaches
the Poisson-Gaussian noise target, after at least two iterations. This discrepancy
principle limits noise fitting; `--no-discrepancy-stop` retains a fixed-iteration
ablation. Spatially varying/non-Gaussian PSFs, distortion, tiled memory use, and
subpixel star registration remain future work. Scale 1 is the memory-conscious
default for joint reconstruction; scale 2 is available explicitly for dithered
datasets.

## Status

This repository is an early but runnable prototype. It includes RAW/DNG
inspection, normalized linear CFA loading, bias/dark/flat calibration against
single master frames or robust in-memory masters built from directories,
CFA-safe star diagnostics, phase-aware drizzle accumulation, several non-neural
demosaic candidates, frequency-guided CFA risk analysis, inverse chroma
refinement, robust joint multi-frame CFA reconstruction, per-channel confidence
maps, TIFF/JPEG export, diagnostic map export, and a minimal Qt preview
GUI with real bias/dark/flat calibration, robust directory masters, and
image/alias-risk/residual/sensor-defect overlays, RAW white-balance selection,
and camera-RGB or linear-sRGB output.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and
[docs/ROADMAP.md](docs/ROADMAP.md).
