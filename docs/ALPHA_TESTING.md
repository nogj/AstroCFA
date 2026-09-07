# AstroCFA Alpha Testing

AstroCFA is ready for technical alpha testing of RAW/DNG loading, CFA-aware
calibration, non-neural demosaicing, TIFF/JPEG export, and diagnostic maps.
The Qt app can now preview a developed frame and switch between image,
alias-risk, and remosaicing-residual overlays.

Useful smoke tests:

```bash
astrocfa-nogui benchmark-debayer

astrocfa-nogui benchmark-debayer --export-prefix benchmark/scene

astrocfa-nogui inspect light.dng --linear-cfa --frequency-cfa

astrocfa-nogui develop light.dng \
  --method inverse-refine \
  -o preview.jpg \
  --preview-stretch astro \
  --export-alias-risk alias-risk.tif \
  --export-residual-map residual.tif

astrocfa-nogui develop light.dng \
  --bias-dir bias/ \
  --dark-dir darks/ \
  --flat-dir flats/ \
  --method inverse-refine \
  -o calibrated.tif

astrocfa-nogui stack light-01.dng light-02.dng light-03.dng \
  --joint-reconstruct \
  --offset 0,0 --offset 0.42,-0.18 --offset -0.31,0.27 \
  -o joint-linear.tif \
  --export-confidence joint-confidence.tif
```

Good reports include:

- operating system and build source;
- `benchmark-debayer` output when possible;
- camera model and RAW format;
- exact command or GUI steps;
- console output with residual/chroma metrics;
- screenshots or crops from the GUI overlay selector;
- whether the alias-risk and residual maps look plausible;
- crops around undersampled stars, saturated stars, and flat-field defects.
- whether `inverse-refine` improves or harms star color compared with
  `inverse-refine-no-star-guard` in `benchmark-debayer`.
- for dithered sequences, initial/final CFA RMSE, robust outlier count, RGB
  direct coverage, and whether the joint confidence map matches the dither pattern.

Do not treat current output as color-managed final data. TIFF export is linear
sensor RGB, and JPEG export is a preview path.
