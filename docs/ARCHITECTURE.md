# AstroCFA Architecture

## Scope

AstroCFA is an offline astrophotography reconstruction engine. It should remain
small, deterministic, and measurable. The project is allowed to be slow when a
slow path produces a more defensible reconstruction.

The primary inputs are RAW camera files, linear DNGs produced by HDRMerge-ng,
and calibrated intermediate frames. The primary outputs are EXR, FITS, TIFF,
and diagnostic maps.

## Pipeline

```text
RAW / DNG CFA frames
  |
  v
linearization
  black subtraction, white level, saturation masks, CFA phase
  |
  v
calibration
  bias, dark, flat, dark-flat, hot/dead pixels, amp glow hints
  |
  v
registration
  star detection, subpixel alignment, optional distortion model
  |
  v
integration
  weighted mean, sigma/winsor clipping, local normalization, CFA drizzle
  |
  v
reconstruction
  baseline demosaic candidates, remosaicing-constrained refinement
  |
  v
astro processing
  background model, PSF estimation, deconvolution, color calibration
  |
  v
output
  scene-linear master, stretched view, uncertainty and residual maps
```

## Core Reconstruction Model

AstroCFA should optimize a scene-linear image against the measured CFA samples:

```text
min X {
  data_fidelity(M(X), Y, noise)
  + luma_regularization(X)
  + conservative_chroma_regularization(X)
  + star_psf_regularization(X)
}
```

Where:

- `Y` is the measured linear CFA image.
- `X` is the reconstructed scene-linear RGB image.
- `M` is the forward model: optics, pixel integration, CFA sampling, noise, and
  quantization where known.
- `data_fidelity` uses a Poisson-Gaussian likelihood rather than plain MSE.

## Noise-Weighted Data Fidelity

AstroCFA's residuals should be interpreted through a sensor noise model. The
initial model is intentionally simple:

```text
variance = read_noise^2 + shot_noise_scale * signal + quantization_noise^2
```

This gives the reconstruction engine a physically meaningful way to decide
whether a remosaicing mismatch is significant. The same numeric difference is
more suspicious in smooth faint sky than in a bright noisy sample.

## Calibration

Calibration operates on linear CFA samples before demosaicing. The current
implementation accepts already-built master bias, dark, and flat frames, or
builds robust in-memory masters from directories of RAW/DNG files. Master inputs
must have the same dimensions and Bayer phase as the light. This is
intentionally strict: a wrong phase or resized master would silently damage
color reconstruction.

The default dark-frame convention is the usual astro one: a master dark already
contains bias, so providing both dark and bias subtracts the dark only. If a
bias-free dark is used, `--dark-excludes-bias` enables explicit bias subtraction
as well.

Flat correction is normalized per CFA phase rather than globally:

```text
calibrated = (light - dark_or_bias) / (flat / mean_flat_for_same_phase)
```

Per-phase flat normalization avoids baking color balance changes into the flat
field step and keeps later white balance/color matrix decisions explicit.

Directory-built masters use a per-pixel robust median with optional MAD clipping
when enough samples are present. This is the right first default for calibration
data because cosmic rays, unstable hot pixels, and occasional bad frames should
not define a master. Persistent CFA master export is deferred until FITS/EXR
support lands, so current directory masters are built in memory for the run.

## Baseline Reconstruction

AstroCFA keeps deliberately explicit demosaic baselines. Bilinear is the lowest
reference. The MHC-like candidate adds Laplacian color-difference correction
while preserving measured CFA samples exactly. The residual-interpolation
candidate reconstructs missing red/blue values by interpolating color
differences against an estimated green plane, which is a better fit for astro
than raw channel interpolation because star structure is mostly luminance while
color should vary more smoothly. These are still stepping stones for RCD,
LMMSE, drizzle, and inverse-refinement work.

## Frequency CFA Diagnostics

AstroCFA includes a frequency-oriented Bayer diagnostic inspired by
luma/chroma multiplexing analysis. The first implementation measures energy in
the CFA carrier patterns `(-1)^x`, `(-1)^y`, and `(-1)^(x+y)` over local tiles.

This is not yet a full Fourier demosaicer. Its purpose is to produce an
alias-risk signal before demosaicing, especially around undersampled stars,
fine diffraction structure, and high-frequency nebulosity where chroma may be
ambiguous.

## Inverse Chroma Refinement

The first inverse-refinement stage keeps the forward model deliberately simple:
the measured CFA channel at each pixel is restored exactly after every
iteration. The free variables are the unmeasured red/blue channels, represented
as chroma residuals against green (`R-G` and `B-G`).

The refinement starts from the frequency-guided demosaic and iteratively updates
interpolated chroma with an edge-aware local solver:

```text
chroma_next = blend(chroma_current, weighted_neighbor_chroma)
```

Weights prefer neighbors across small luminance changes, give extra authority
to real measured red/blue anchors, downweight invalid or clipped samples, and
increase chroma caution in high alias-risk tiles.

This is a conservative first inverse stage rather than a full variational
optimizer. Its purpose is to make the underdetermined part of demosaicing
explicit and measurable: preserve sensor samples, reduce unsupported chroma
roughness, and prepare the codebase for later Poisson-Gaussian likelihood,
PSF-aware star terms, and CFA-drizzle-informed confidence maps.

## Internal CFA Representation

The internal CFA frame stores linear normalized sample values plus validity and
clipping flags. This matters for astro reconstruction because saturated stellar
cores, invalid registration borders, and rejected pixels should not contribute
to the same data-fidelity term as trustworthy sensor samples.

```text
CfaSample = value + valid + clipped
```

The first loader target subtracts per-phase black levels and normalizes by the
camera white level without applying demosaicing, white balance, gamma, or color
conversion.

## CFA-Safe Star Detection

Early star diagnostics operate on a 2x2 CFA-safe luminance proxy. Each proxy
sample is built from measured linear CFA values in one Bayer cell, skipping
invalid samples and excluding clipped cells from background estimation. This
keeps star detection anchored to sensor measurements before RGB interpolation.

The initial detector estimates sky background with median and MAD rather than
ordinary mean and standard deviation. That is a better default for sparse stars,
hot pixels, weak nebulosity, and frames where bright structures should not
define the detection threshold.

## Star Registration

Registration belongs to AstroCFA only when it is star-based. HDRMerge-ng owns
general exposure-bracket alignment; AstroCFA owns light-frame alignment driven
by stellar signal.

The initial implementation estimates integer translations on the same CFA-safe
2x2 luminance proxy used for star diagnostics. It is intentionally a coarse
stage and uses a sparse sampling grid by default so full-resolution RAWs remain
interactive enough for inspection. Later stages should refine this with
centroids, matched star lists, subpixel offsets, rotation, optical distortion,
and per-frame quality weights.

## CFA Drizzle

When lights are dithered, AstroCFA should avoid demosaicing each exposure first.
Instead, it should register calibrated CFA frames and accumulate same-phase
samples into a higher-information reconstruction grid. This allows real samples
from different frames to reduce Bayer ambiguity before interpolation.

```text
calibrated Bayer frames + subpixel offsets
  -> per-phase drizzle accumulation
  -> inverse RGB reconstruction
```

This is the first major differentiator from ordinary RAW development and simple
astro scripts: the mosaic is treated as sensor data until the stack has gathered
as much information as the capture session can provide.

The first accumulator stores separate high-resolution planes for `R`, `G1`,
`B`, and `G2`. Valid, unclipped samples are deposited with subpixel offsets and
frame weights. Coverage and accumulated weight are first-class outputs because
they tell the reconstruction stage which colors/frequencies are measured and
which remain inferred.

## Diagnostic Maps

Every high-quality run should be able to emit:

- clipping mask;
- hot/dead pixel mask;
- noise sigma map;
- remosaicing residual;
- chroma confidence map;
- aliasing/undersampling map;
- star saturation and halo risk map;
- integration rejection map.

These maps are not secondary polish. They are part of the product identity:
AstroCFA should show where the result is measured and where it is inferred.

## Positioning Against Existing Astro Pipelines

Existing astro tools already cover calibration, registration, stacking, debayer,
stretching, and color work. AstroCFA should not try to win by having more
general-purpose features.

The differentiator is stricter:

```text
Every reconstructed RGB pixel should remain accountable to measured CFA samples.
```

That means AstroCFA should expose diagnostics that most ordinary user-facing
pipelines hide:

- per-phase clipping before demosaicing;
- per-phase flat-field sanity checks;
- G1/G2 split for CFA balance and flat calibration diagnostics;
- same-phase drizzle accumulation coverage;
- RGB-to-CFA remosaicing residual;
- chroma ambiguity around undersampled stars;
- PSF change before and after reconstruction/deconvolution.

The desired behavior is not always the sharpest-looking output. In undersampled
stellar fields, a faithful result may suppress unsupported chroma or flag
uncertainty instead of producing attractive but unmeasured color detail.

## Non-Goals

- Catalog management.
- Preset-based creative editing.
- Neural denoising, neural demosaicing, or generative detail recovery.
- In-place replacement for Lightroom, darktable, RawTherapee, Siril, or
  PixInsight.
- A large GUI before the reconstruction engine is validated.
