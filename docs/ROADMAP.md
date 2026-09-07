# AstroCFA Roadmap

## Milestone 0: Project Skeleton

- CMake project.
- Qt GUI entry point as `astrocfa`.
- CLI entry point as `astrocfa-nogui`.
- Architecture and product thesis.
- Basic tests for help/version commands.

## Milestone 1: RAW Inspection

- Read RAW/DNG metadata through LibRaw.
- Report dimensions, active area, black/white levels, CFA pattern, ISO,
  exposure time, camera model, and clipping statistics.
- Report saturated component topology before demosaicing for star-core handling.
- Preserve orientation and CFA phase information.
- Add synthetic fixtures for Bayer phase handling.

## Milestone 2: Linear CFA Frame

- Internal `CfaFrame` type with per-sample value and phase.
- Black subtraction and white normalization.
- Basic FITS/EXR export for linear CFA diagnostics.
- Unit tests for RGGB/BGGR/GRBG/GBRG addressing.
- Remosaicing residual for validating reconstructed RGB against measured CFA.

## Milestone 3: Calibration Frames

- Bias, dark, flat, and dark-flat application.
- Hot/dead pixel detection from calibration masters, with same-phase repair and
  an exportable defect map. Implemented.
- Robust master frame creation.
- Rejection maps for cosmic rays and unstable pixels.

## Milestone 4: Registration And Integration

- Star detection on CFA-safe luminance proxy.
- Subpixel translation registration.
- Weighted stacking with sigma or winsorized clipping.
- CFA-aware drizzle accumulation.
- Robust joint CFA reconstruction from translated lights, with per-channel
  confidence, transient rejection, and known Gaussian PSF per light. Implemented.
- Deterministic multi-frame ablation benchmark with known dithers, sensor noise,
  transients, photometry, FWHM, elongation, and RGB/chroma truth metrics. Implemented.
- Robust per-frame Gaussian PSF estimation directly from CFA star samples and
  relative-seeing forward models. Implemented.
- Select iterations by the Poisson-Gaussian discrepancy principle. Implemented.
- Extend the retained per-star measurements into spatially varying,
  non-Gaussian PSF fields.

## Milestone 5: Baseline Reconstruction

- RCD-style demosaic baseline for Bayer.
- LMMSE-style noisy-frame mode.
- Remosaicing residual metric.
- Synthetic star-field test images with known ground truth.

## Milestone 6: Faithful Astro Reconstruction

- Tile-based inverse refinement initialized from baseline demosaic.
- Poisson-Gaussian data fidelity.
- Conservative chroma regularization.
- Star-aware PSF preservation.
- Chroma confidence and aliasing maps.

## Milestone 7: Astro Development

- Source-resistant, robust polynomial background model with retained sky level,
  optional neutralization, and diagnostic surface export. Implemented.
- User masks and multiscale protection for fields dominated by extended signal.
- Spatially varying/non-Gaussian PSF estimation.
- Regularized Richardson-Lucy and Wiener deconvolution.
- RAW as-shot/daylight/custom white balance and camera-to-linear-sRGB matrix
  conversion. Implemented.
- ICC/DCP camera profiles and wide-gamut linear working spaces.
- Photometric color calibration hooks.
- Luminance-ratio arcsinh and published generalized hyperbolic stretch with
  exposure, levels, color preservation, and highlight protection. Implemented.
- Interactive histogram and reusable development presets.
