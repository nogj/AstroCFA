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
- Source-resistant additive background modeling with an auditable surface map.
- Luminance-ratio arcsinh and Generalized Hyperbolic Stretch with star highlights.
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
selecting an input, loading bias/dark/flat masters, previewing the reconstruction,
switching diagnostic overlays, and exporting
TIFF/JPEG output. A calibration selector accepts either one master RAW or a
directory of RAW frames; directories are combined into a robust in-memory master.
The reconstruction is purpose-specific: AstroCFA separates diffuse structure and
point sources directly on the CFA, learns a shared stellar ePSF, and selects a
radial, achromatic 2D, or chromatic 2D model by held-out CFA evidence. White balance
and output-space selectors keep sensor RGB, as-shot/daylight metadata, and linear
sRGB conversion explicit.

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
  -o calibrated-preview.jpg --preview-stretch astro
astrocfa-nogui calibrate light.dng --dark-dir darks/ --flat-dir flats/ \
  -o calibrated-preview.jpg --preview-stretch astro
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
astrocfa-nogui develop input.dng -o preview.jpg \
  --preview-stretch astro --export-alias-risk alias.tif --export-residual-map residual.tif
astrocfa-nogui develop input.dng --white-balance daylight \
  --output-space srgb -o developed-linear-srgb.tif
astrocfa-nogui develop input.dng --background gradient \
  --export-background background.tif -o gradient-corrected.tif
astrocfa-nogui develop input.dng --tone ghs \
  --stretch-factor 3.0 --local-intensity 8 --symmetry-point 0.08 \
  --protect-highlights 0.80 --saturation 1.05 -o developed-ghs.tif
astrocfa-nogui develop input.dng \
  --wb-multipliers 2.1,1.0,1.4 --output-space camera -o developed-camera-rgb.tif
astrocfa-nogui develop light.dng --bias master-bias.dng --dark master-dark.dng \
  --flat master-flat.dng -o calibrated.tif \
  --export-defect-map sensor-defects.tif
```

`develop` writes scene-linear 16-bit RGB TIFF output or an 8-bit JPEG preview.
Its default color mode uses as-shot WB and the camera-to-sRGB matrix when the RAW
provides them, with explicit daylight, unity, custom-multiplier, and camera-RGB
alternatives. JPEG output can use an astro-oriented arcsinh preview stretch while
TIFF remains linear. `develop` and `calibrate` always use AstroCFA's adaptive
morphological reconstruction. The CLI and GUI report the selected stellar model,
held-out scores, source counts, and chromatic ePSF transform, so model complexity
is visible rather than hidden behind a preset. Established demosaicers remain only
as internal benchmark baselines.

Background correction is deliberately opt-in. `--background gradient` fits an
additive degree-2 surface from low-luminance samples in spatial tiles, then uses
asymmetric robust weighting to suppress stars, halos, and positive extended
signal. It removes only the modeled spatial variation and preserves the median
sky level per channel. `--background neutral` additionally makes that retained
level achromatic; `--export-background` exposes the fitted surface for review.

Tone mapping is also opt-in in the CLI. `--tone linear|arcsinh|ghs` applies
exposure, levels, and a selected stretch to TIFF or JPEG output. Arcsinh and GHS
operate on linear luminance and rescale RGB by the luminance ratio, preserving
stellar color until gamut compression is required. GHS exposes its local
intensity, symmetry point, and linear shadow/highlight protection segments.
Reports include selected black/white points, clipping, and per-pixel gamut
compression. Without `--tone`, TIFF remains scene-linear; the older
`--preview-stretch astro` remains a JPEG-only shortcut.

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

The experimental `cfa-mca-poc` row performs direct, noise-weighted CFA fitting
with separate diffuse and point-source components. It can learn a shared radial
PSF from a subset of stars selected for profile-bin information gain, while
`cfa-mca-independent-poc` is the no-sharing ablation. Both remain benchmark-only
research paths, not user-facing development methods. Controlled fixtures can use
`--stars`, `--common-star-sigma`, `--moffat-beta`, `--psf-ellipticity`,
`--psf-angle`, `--chromatic-psf-shift`, `--chromatic-psf-scale`,
`--star-flux-scale`, and `--star-phase` to test the mechanism without changing
the historical fixture. A flux-normalized, oversampled 2D ePSF is also fitted by
alternating closed-form RGB flux estimates with constrained profile updates;
source positions are then reoptimized against the ePSF and the diffuse component
is conservatively re-estimated after subtracting the fitted stars. A low-rank
chromatic extension can estimate global R/B shifts and scale relative to G,
guarded by both CFA holdout and a six-parameter BIC penalty.
`cfa-mca-no-epsf-poc` and `cfa-mca-achromatic-epsf-poc` are explicit ablations.
The shared row also reports the selected model, source counts, full-matrix
log-determinant information gain, and normalized CFA holdout scores for all
five candidate families. A second morphology-aware detection pass is only used
when unseen CFA samples select the ePSF.

The benchmark can also export a standards-based linear CFA DNG fixture and score
external RGB TIFF candidates through repeated `--external name=path` arguments;
`--external-srgb` performs a declared standard sRGB transfer decode first.
`--csv` writes the internal and external rows in one machine-readable
table. The exported manifest and star catalogue record the exact test contract;
see [benchmarks/EXTERNAL_DEBAYER.md](benchmarks/EXTERNAL_DEBAYER.md) for the
RawTherapee and Siril comparison workflow.
The first RawTherapee 5.13 baseline is reported candidly in
[benchmarks/RESULTS.md](benchmarks/RESULTS.md).

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
camera-RGB or linear-sRGB output, and an optional protected background model.
The GUI tone controls export the same arcsinh, GHS, or linear-level rendering
shown in its preview.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and
[docs/ROADMAP.md](docs/ROADMAP.md).
