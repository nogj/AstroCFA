# AstroCFA Alpha Testing

AstroCFA is ready for technical alpha testing of RAW/DNG loading, CFA-aware
calibration, non-neural demosaicing, TIFF/JPEG export, and diagnostic maps.

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
```

Good reports include:

- operating system and build source;
- `benchmark-debayer` output when possible;
- camera model and RAW format;
- exact command or GUI steps;
- console output with residual/chroma metrics;
- whether the alias-risk and residual maps look plausible;
- crops around undersampled stars, saturated stars, and flat-field defects.

Do not treat current output as color-managed final data. TIFF export is linear
sensor RGB, and JPEG export is a preview path.
