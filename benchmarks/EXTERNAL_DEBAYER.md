# External Astro Debayer Benchmark

AstroCFA can score an external demosaicer with the same reconstruction metrics
used for its internal candidates. The benchmark is deliberately linear and
camera-space: it measures reconstruction, not rendering preferences.

## 1. Export A Fixture

```bash
astrocfa-nogui benchmark-debayer \
  --width 512 --height 384 --seed 7 --noise astro \
  --export-prefix benchmark/scene-007
```

This writes:

- `scene-007-input.dng`: 16-bit linear Bayer input, RGGB by default.
- `scene-007-truth.tif`: 16-bit linear RGB ground truth.
- `scene-007-cfa.tif`: grayscale visualization of the mosaic.
- `scene-007-stars.csv`: exact synthetic star catalogue.
- `scene-007-manifest.json`: fixture parameters and candidate contract.
- JPEG previews for AstroCFA's internal methods.

The DNG contains unity as-shot white balance, explicit CFA metadata, zero black,
and a 65535 white level. It is intentionally uncompressed for broad RAW decoder
compatibility.

## 2. Run An External Demosaicer

The candidate output must be a top-left, full-size RGB TIFF with unsigned 8/16-bit
or float32 samples. Values must be finite and in `[0, 1]`. For a valid comparison:

- keep the image scene-linear;
- use unity white balance;
- keep camera RGB and disable the input color transform;
- disable exposure compensation, tone curves, denoising and sharpening;
- disable lens correction, cropping and resizing.

RawTherapee can be automated with a versioned neutral `.pp3` profile. In that
profile select the demosaicing method under Bayer Sensor, use unity white balance,
select `No profile` as the input profile, choose a linear-gamma output, and disable
the remaining processing modules. Save separate profiles for RCD, AMaZE, LMMSE,
and any other candidate, then run for example:

```bash
rawtherapee-cli -o benchmark/rt-rcd.tif \
  -p benchmarks/profiles/rawtherapee-rcd-linear.pp3 \
  -b16 -tz -Y -c benchmark/scene-007-input.dng
```

Processing-profile keys change between RawTherapee releases, so commit the actual
profile together with its `rawtherapee-cli --version` output. Do not treat a
profile made for an untested release as equivalent.

Siril can import and debayer the DNG through `convertraw ... -debayer`. Export the
result as a linear TIFF without color calibration or stretching before scoring.
Record the Siril version, debayer method and conversion script next to the result.

## 3. Score Every Candidate

Candidates are repeatable `name=path` arguments. The scene options must match the
fixture manifest:

```bash
astrocfa-nogui benchmark-debayer \
  --width 512 --height 384 --seed 7 --noise astro \
  --external rawtherapee-rcd=benchmark/rt-rcd.tif \
  --external rawtherapee-amaze=benchmark/rt-amaze.tif \
  --external siril=benchmark/siril.tif \
  --csv benchmark/results.csv
```

Lower values are better. No single metric is sufficient: inspect global RGB and
chroma errors together with stellar false color, luminance, aperture flux, FWHM,
elongation, and CFA residual. A method that sharpens stars while biasing flux is
not an unconditional improvement.

## Reporting Rules

- Publish the fixture manifest, external profiles/scripts and software versions.
- Run multiple seeds in both `none` and `astro` noise modes.
- Do not tune a method on the same seeds used for the final comparison.
- Report every metric, including regressions and CFA residual.
- Keep visual examples secondary to the numeric table.
