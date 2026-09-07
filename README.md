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
selecting an input, choosing a reconstruction mode, and launching the major
astro stages.

`astrocfa-nogui` is the command-line executable for scripts, batch processing,
and reproducible long-running workflows.

## Initial CLI Shape

```bash
astrocfa-nogui inspect input.dng
astrocfa-nogui inspect input.dng --linear-cfa --star-candidates --noise-model --frequency-cfa
astrocfa-nogui benchmark-debayer
astrocfa-nogui benchmark-debayer --width 192 --height 128 --seed 7 --export-prefix bench/astro
astrocfa-nogui calibrate light.dng --dark master-dark.dng --flat master-flat.dng \
  --method inverse-refine -o calibrated-preview.jpg --preview-stretch astro
astrocfa-nogui calibrate light.dng --dark-dir darks/ --flat-dir flats/ \
  --method inverse-refine -o calibrated-preview.jpg --preview-stretch astro
astrocfa-nogui stack light1.dng light2.dng light3.dng --cfa-drizzle --scale 2 \
  --offset 0,0 --offset 0.42,-0.18 --offset -0.31,0.27
astrocfa-nogui stack light1.dng light2.dng light3.dng --cfa-drizzle --auto-register
astrocfa-nogui develop input.dng --method bilinear-baseline -o baseline.tif
astrocfa-nogui develop input.dng --method malvar-baseline -o malvar.tif
astrocfa-nogui develop input.dng --method residual-interpolation -o ri.tif
astrocfa-nogui develop input.dng --method frequency-guided -o preview.jpg
astrocfa-nogui develop input.dng --method inverse-refine -o preview.jpg \
  --preview-stretch astro --export-alias-risk alias.tif --export-residual-map residual.tif
astrocfa-nogui develop light.dng --bias master-bias.dng --dark master-dark.dng \
  --flat master-flat.dng --method inverse-refine -o calibrated.tif
```

`develop` writes scene-linear 16-bit RGB TIFF output or an 8-bit JPEG preview.
JPEG output can use an astro-oriented arcsinh preview stretch while TIFF remains
linear. The `inverse-refine` method starts from frequency-guided residual
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
error, and remosaicing residual. This is the early guardrail for keeping
AstroCFA's reconstruction work measurable rather than merely aesthetic.

## Status

This repository is an early but runnable prototype. It includes RAW/DNG
inspection, normalized linear CFA loading, bias/dark/flat calibration against
single master frames or robust in-memory masters built from directories,
CFA-safe star diagnostics, phase-aware drizzle accumulation, several non-neural
demosaic candidates, frequency-guided CFA risk analysis, inverse chroma
refinement, TIFF/JPEG export, and diagnostic map export.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and
[docs/ROADMAP.md](docs/ROADMAP.md).
