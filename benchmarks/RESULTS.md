# External Debayer Results

## RawTherapee 5.13 Baseline

This is an initial engineering baseline, not a state-of-the-art claim. It uses
one deterministic 512x384 synthetic scene (`seed=7`) with five stars. Every
RawTherapee output was generated as a 16-bit sRGB TIFF using the checked-in 5.13
profiles and decoded with AstroCFA's `--external-srgb` path before measurement.
The complete machine-readable tables are
[`rawtherapee-5.13-seed7-noiseless.csv`](results/rawtherapee-5.13-seed7-noiseless.csv)
and [`rawtherapee-5.13-seed7-astro.csv`](results/rawtherapee-5.13-seed7-astro.csv).

Selected metrics are shown below. Lower is better.

### Noiseless CFA

| Method | RGB RMSE | Chroma MAE | Star false color | Star luma RMSE | Flux relative | FWHM relative |
|---|---:|---:|---:|---:|---:|---:|
| AstroCFA inverse-refine | 0.000845 | 0.000130 | 0.020722 | 0.036036 | 0.135720 | 0.151752 |
| RawTherapee RCD | 0.001020 | 0.000014 | 0.017044 | 0.049078 | 0.124122 | 0.068324 |
| RawTherapee AMaZE | 0.001090 | 0.000020 | 0.025119 | 0.044029 | 0.080966 | 0.038909 |
| RawTherapee LMMSE | **0.000439** | **0.000012** | **0.008606** | **0.016050** | **0.027406** | **0.016970** |

### Poisson-Gaussian Noise And Hot Pixels

| Method | RGB RMSE | Chroma MAE | Star false color | Star luma RMSE | Flux relative | FWHM relative |
|---|---:|---:|---:|---:|---:|---:|
| AstroCFA inverse-refine | 0.037699 | 0.003866 | 0.020871 | 0.036557 | 0.131817 | 0.161631 |
| RawTherapee RCD | **0.035482** | 0.004688 | 0.019496 | 0.045111 | 0.113988 | 0.085993 |
| RawTherapee AMaZE | 0.038301 | 0.004495 | 0.021580 | 0.028165 | 0.041032 | 0.083727 |
| RawTherapee LMMSE | 0.037594 | **0.002738** | **0.008992** | **0.013213** | **0.018649** | **0.083551** |

LMMSE is the strongest method on most stellar metrics in both variants. In the
noisy case, RCD has the lowest RGB RMSE among the methods in this table, while
the simple AstroCFA bilinear baseline is lower still (`0.033577`) because global
error is dominated by the smooth noisy background. This exposes a benchmark
design issue: global RMSE alone rewards smoothing and must not rank astro
reconstruction by itself.

AstroCFA's current advantage is exact measured-sample fidelity. Its current
inverse refinement is not competitive enough in stellar photometry or FWHM on
this fixture. The next algorithmic work should target that measured gap and the
next benchmark revision must include more stars, multiple PSF widths, multiple
seeds, and signal-stratified background/source metrics.
