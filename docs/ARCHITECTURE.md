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

Cosmetic correction is derived exclusively from calibration masters. A dark
sample is considered hot when it exceeds its same-phase local median by both a
robust sigma threshold and an absolute signal floor. A flat sample is considered
dead when its local same-phase response ratio and absolute deficit are both low.
This dual criterion limits false positives across amp glow and flat-field
gradients. Defective samples are excluded from flat phase normalization and are
replaced after calibration by a median of valid, non-defective neighbors from
the same CFA phase. If too few such measurements exist, the sample is marked
invalid instead of being fabricated.

The resulting defect map records hot, dead, and invalid-master evidence
separately. It can be inspected in the GUI or exported by the CLI, and the run
report records detected, repaired, and unrepaired counts.

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

The solver also includes a star chroma guard. Compact luminance peaks are where
single-frame Bayer demosaicing most often invents magenta/green stellar cores.
For unmeasured red/blue channels, AstroCFA can reduce unsupported chroma in
these PSF-like peaks while still restoring the measured CFA channel exactly.
This is intentionally tunable: white undersampled stars benefit from stronger
guarding, while fields with reliable color sampling or larger PSFs can use a
lighter setting.

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

## Joint Multi-Frame CFA Reconstruction

The first joint solver estimates one RGB image from all calibrated lights rather
than demosaicing and stacking each light independently. Its forward model samples
the observed color channel from the common output grid at each frame's translated
sensor coordinate. The residual is weighted by the Poisson-Gaussian noise model
and a Huber influence function, then backprojected through the same bilinear
kernel. This provides robust resistance to transient samples while retaining an
auditable residual against every input measurement.

Initialization combines phase-separated measurements on the common grid and
uses the first frame's conservative reconstruction only to fill unconstrained
locations. Chroma and green-luminance regularization are edge-aware. Direct CFA
support is immutable in the single-frame case. With multiple noisy lights it is
a weighted fidelity constraint rather than a hard lock, allowing the common
latent estimate to denoise conflicting observations. The solver exports
normalized R/G/B support as a confidence image and reports coverage, outlier
count, normalized error, initial/final CFA RMSE, and robust reduced chi-square.
By default, iteration stops once that chi-square is consistent with the declared
Poisson-Gaussian noise model. A minimum iteration count prevents the initialized
estimate from terminating the solve before robust rejection has begun.

The deterministic `benchmark-joint` harness compares demosaic-then-average,
phase-aware initialization, PSF-aware least-squares joint inversion, Huber-robust
joint inversion without PSF information, auto-estimated PSF inversion, and
PSF-aware robust inversion against known RGB truth. It reports global RGB/chroma
errors plus stellar aperture flux, FWHM, elongation, false color, and luminance errors. Its
variable-seeing mode is an explicit ablation of the PSF operator and estimator.

Each light may provide a Gaussian PSF sigma in sensor pixels. The forward model
combines that convolution with subpixel bilinear sampling in one normalized
stencil; residual backprojection uses its exact transpose. PSF updates use a
conservative SIRT normalization, flux-conserving bilateral luminance diffusion,
and the existing chroma prior. Translation-only inputs retain the earlier
diagonal update and exact single-frame CFA behavior.

Automatic PSF estimation detects isolated, unsaturated stars on the CFA-safe
luminance proxy, rejects low-SNR and elongated candidates, then fits a circular
Gaussian directly to original Bayer samples. Separate per-channel amplitude and
background nuisance terms avoid treating stellar color as profile structure.
Robust median/MAD aggregation produces one seeing estimate per frame. Since the
latent image already contains the sharpest observation's optical blur, only the
additional variance relative to the sharpest frame is passed to the solver, with
a conservative default strength of 0.75. Individual star measurements are kept
for auditability and a future spatial PSF field. The model does not yet cover
non-Gaussian/spatially varying PSFs, lens distortion, photometric scale,
background offsets, tiled execution, or inverse-Hessian uncertainty.

## Protected Background Modeling

Astro background correction operates on reconstructed linear camera RGB before
white balance and color conversion. The image is divided into tiles and only a
low-luminance quantile from each tile contributes to a representative RGB
sample. This lower-envelope sampling reduces direct contamination from stars
without blurring or resampling the science image. If a small image would provide
too few samples, the tile size is reduced automatically until the requested
polynomial is constrained.

A degree-0, degree-1, or degree-2 spatial surface is fit independently to each
camera channel while sharing robust weights derived from luminance residuals.
Positive residuals are rejected more aggressively than negative residuals:
stars, halos, and nebulosity should not pull the instrumental background upward.
The fitted surface is subtracted additively and its median per-channel level is
restored, preserving the sky pedestal unless neutralization is explicitly
selected.

The model reports tile counts, downweighted samples, robust residual scale,
preserved RGB levels, and per-channel peak-to-peak gradients. It also returns
the complete fitted surface as a diagnostic image. Correction is opt-in because
no source mask can prove that a smooth structure is instrumental rather than
astronomical. Future work should add user masks, multiscale segmentation, and
comparison against dither-derived sky models.

## RAW Color Development

Demosaicing produces linear camera RGB. Color development is a separate stage
so white balance and output-space choices cannot alter remosaicing residuals or
the reconstruction solver. AstroCFA reads LibRaw's as-shot and daylight channel
multipliers plus its camera-to-sRGB matrix. White-balance multipliers are
normalized around green and applied before the matrix.

Automatic mode prefers as-shot metadata, then daylight metadata, then sensor
unity. It applies the camera matrix only when a finite nonzero matrix is present.
Explicit requests fail when their required metadata is absent instead of silently
substituting a different color interpretation. A custom RGB multiplier path and
camera-RGB output preserve reproducible alternatives for narrowband, modified
cameras, and external profiling workflows.

The transform remains scene-linear and does not clip internally. Negative and
over-range pixels are counted in the development report. TIFF/JPEG encoding
currently clamps to the integer output range; the arcsinh stretch and display
gamma belong only to the JPEG/GUI preview path. ICC embedding, camera-profile
interpolation, working spaces wider than sRGB, and photometric color calibration
remain future work.

## Astro Tone Development

The tone engine supports linear levels, arcsinh, and the published Generalized
Hyperbolic Stretch family. GHS uses `D = exp(stretch_factor) - 1`, local intensity
`b`, a symmetry point, and tangent-linear shadow/highlight protection segments,
then normalizes the complete function to map zero and one exactly. The equations
follow the [GHS process reference](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html#transformation_equations).

Nonlinear transforms operate in a color-preserving mode. The engine derives
linear luminance, transforms it once, and scales RGB by the transformed/original
luminance ratio. Optional saturation is applied around the transformed luminance.
If a channel would leave the display gamut, chroma is compressed toward that
luminance rather than independently clipping the channel. This is especially
important for maintaining star color through aggressive stretches.

Automatic black and white points use configurable luminance percentiles. Manual
points, scene-linear exposure compensation, and clipping counts remain explicit.
The report records shadow/highlight clipping and gamut-compressed pixels. CLI
tone mapping is opt-in so an ordinary TIFF stays linear; the Qt preview and its
export deliberately share the selected tone transform.

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
