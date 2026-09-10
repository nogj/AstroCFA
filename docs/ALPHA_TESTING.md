# AstroCFA Alpha Testing

AstroCFA is ready for technical alpha testing of RAW/DNG loading, CFA-aware
calibration, non-neural demosaicing, TIFF/JPEG export, and diagnostic maps.
The Qt app can now preview a developed frame and switch between image,
alias-risk, and remosaicing-residual overlays.

Useful smoke tests:

```bash
astrocfa-nogui benchmark-debayer

astrocfa-nogui benchmark-debayer --export-prefix benchmark/scene
astrocfa-nogui benchmark-debayer --width 512 --height 384 --seed 7 \
  --export-prefix benchmark/scene-007
astrocfa-nogui benchmark-debayer --width 512 --height 384 --seed 7 \
  --external-srgb rawtherapee-rcd=benchmark/rt-rcd.tif \
  --csv benchmark/results.csv

astrocfa-nogui benchmark-joint --export-prefix benchmark/joint

astrocfa-nogui benchmark-joint --seeing variable --transients 4

astrocfa-nogui inspect light.dng --linear-cfa --frequency-cfa --estimate-psf

astrocfa-nogui develop light.dng \
  -o preview.jpg \
  --preview-stretch astro \
  --export-alias-risk alias-risk.tif \
  --export-residual-map residual.tif

astrocfa-nogui develop light.dng \
  --white-balance daylight \
  --output-space srgb \
  -o daylight-linear-srgb.tif

astrocfa-nogui develop light.dng \
  --background gradient \
  --export-background fitted-background.tif \
  -o gradient-corrected.tif

astrocfa-nogui develop light.dng \
  --tone ghs --stretch-factor 3.0 --local-intensity 8 \
  --symmetry-point 0.08 --protect-highlights 0.80 \
  -o ghs-render.tif

astrocfa-nogui develop light.dng \
  --bias-dir bias/ \
  --dark-dir darks/ \
  --flat-dir flats/ \
  -o calibrated.tif

astrocfa-nogui stack light-01.dng light-02.dng light-03.dng \
  --joint-reconstruct \
  --offset 0,0 --offset 0.42,-0.18 --offset -0.31,0.27 \
  --psf-sigma 0.72 --psf-sigma 0.91 --psf-sigma 0.68 \
  -o joint-linear.tif \
  --export-confidence joint-confidence.tif

astrocfa-nogui stack light-01.dng light-02.dng light-03.dng \
  --joint-reconstruct --auto-register --auto-psf \
  -o joint-auto-psf.tif \
  --export-confidence joint-auto-confidence.tif
```

Good reports include:

- operating system and build source;
- `benchmark-debayer` output when possible;
- external demosaicer version, processing profile and fixture manifest when using
  `--external`;
- `benchmark-joint` output for dithered-data or joint-solver reports;
- camera model and RAW format;
- exact command or GUI steps;
- console output with residual/chroma metrics;
- screenshots or crops from the GUI overlay selector;
- whether the alias-risk and residual maps look plausible;
- crops around undersampled stars, saturated stars, and flat-field defects.
- the selected AstroCFA point-source model, holdout scores, source counts, and any
  reported red/blue ePSF shift or scale.
- selected WB multipliers, output color space, and negative/out-of-range counts.
- background tile/rejection counts, gradient amplitudes, and the exported fitted
  surface; inspect it for leaked nebulosity or halos before trusting correction.
- tone black/white points, clipping and gamut-compression counts; compare color
  and profile width in bright stars against the linear source.
- for dithered sequences, initial/final CFA RMSE, robust outlier count, RGB
  direct coverage, and whether the joint confidence map matches the dither pattern.
- for variable-seeing sequences, estimated FWHM/star counts per light, relative
  PSF sigmas, reduced chi-square, and whether discrepancy stopping triggered.

Do not treat current output as color-managed final data. TIFF export is linear
sensor RGB, and JPEG export is a preview path.
