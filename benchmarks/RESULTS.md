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

## Phase-Diverse Non-Local Proof Of Concept

This experiment tests a narrower hypothesis: whether repeated stars can provide
enough complementary CFA samples to learn a shared, non-Gaussian PSF without a
neural prior. It uses 25 noisy stars with independent RGB amplitudes and a common
Moffat PSF (`sigma=0.65`, `beta=2.5`). Candidate locations come from unsaturated
CFA blocks, not the truth catalogue. A deterministic 20% CFA holdout chooses
between independent Gaussian fits, a shared Gaussian, and an empirical radial
profile. The profile-learning subset is selected greedily by full-matrix
log-determinant information gain over the coupled radial-bin operator.

```bash
astrocfa-nogui benchmark-debayer --width 512 --height 384 --seed 7 \
  --noise astro --stars 25 --common-star-sigma 0.65 --moffat-beta 2.5 \
  --star-flux-scale 0.5 --star-phase diverse
```

The complete table is
[`rawtherapee-5.13-phase-diverse-moffat.csv`](results/rawtherapee-5.13-phase-diverse-moffat.csv).

| Method | RGB RMSE | Star false color | Star luma RMSE | Flux relative | FWHM relative | Elongation |
|---|---:|---:|---:|---:|---:|---:|
| CFA-MCA independent | 0.001820 | 0.010141 | 0.036772 | 0.161493 | 0.134716 | 0.011822 |
| CFA-MCA radial profile | 0.001724 | **0.007823** | 0.034803 | 0.145904 | 0.093717 | **0.010589** |
| CFA-MCA selected ePSF | **0.001145** | 0.009115 | 0.015820 | **0.071768** | **0.009469** | 0.011445 |
| RawTherapee AMaZE | 0.038298 | 0.013028 | 0.017587 | 0.119721 | 0.076342 | 0.186540 |
| RawTherapee LMMSE | 0.037597 | 0.010345 | **0.014154** | **0.111037** | **0.073422** | 0.176187 |

The radial profile improves every listed stellar metric over its independent
ablation. Subpixel position refinement lets the 2D ePSF improve flux and shape
again: it beats LMMSE and AMaZE in aperture flux, FWHM, elongation, RGB RMSE, and
chroma MAE, and narrowly beats LMMSE in false color. LMMSE retains the best
stellar luminance (`0.014154` versus `0.015820`). Global RMSE remains dominated
by background treatment and must not be interpreted alone as an overall ranking.

The auditable selector detected 64 candidates, validated 23 against unseen CFA
samples, used 12 for radial-profile learning, and reconstructed all 25 injected
stars after the morphology-aware pass. The ePSF won the normalized holdout score
(`1.421883`) over the radial profile (`1.423399`), shared Gaussian (`1.427574`),
and independent fits (`1.429777`). The radial selector's full information gain
over the regularizing prior was `506.346755` in natural-log determinant units.

A matched fixture in which every star center is forced to the same Bayer and
subpixel phase had lower information gain (`493.282205`) but did not perform
worse consistently. Therefore this PoC supports
**shared non-Gaussian PSF estimation**, but does not yet establish the stronger
claim that the current full-matrix selector converts CFA phase
diversity. Full operator conditioning, multiple seeds, and spatially varying
PSFs remain required before this can graduate from a research candidate.

## Flux-Normalized 2D ePSF Proof Of Concept

The next experiment replaces the circular profile with a positive,
flux-normalized 2D effective PSF sampled at twice the sensor resolution. Its RGB
source fluxes are solved in closed form for every profile update. The empirical
surface can therefore represent anisotropy without assigning a separate color
shape to each star. A morphology-aware second pass revisits candidates rejected
by the circular Gaussian seed model, but it runs only after unseen CFA samples
select the ePSF family.

The fixture uses the same 25-star noisy Moffat scene with ellipticity `0.35` and
angle `0.55` radians:

```bash
astrocfa-nogui benchmark-debayer --width 512 --height 384 --seed 7 \
  --noise astro --stars 25 --common-star-sigma 0.65 --moffat-beta 2.5 \
  --psf-ellipticity 0.35 --psf-angle 0.55 --star-flux-scale 0.5 \
  --star-phase diverse
```

The complete internal table is
[`astrocfa-epsf-anisotropic.csv`](results/astrocfa-epsf-anisotropic.csv).

| Method | RGB RMSE | Star false color | Star luma RMSE | Flux relative | FWHM relative | Elongation |
|---|---:|---:|---:|---:|---:|---:|
| Frequency-guided | 0.037047 | 0.018719 | 0.021107 | 0.162777 | 0.058359 | 0.168286 |
| CFA-MCA without ePSF | 0.002902 | 0.019552 | 0.061336 | 0.444819 | 0.311622 | 0.723429 |
| CFA-MCA selected ePSF | **0.001246** | 0.014978 | 0.018463 | **0.081504** | **0.010301** | **0.038724** |
| RawTherapee AMaZE | 0.038299 | 0.013414 | 0.017195 | 0.127884 | 0.064841 | 0.184590 |
| RawTherapee LMMSE | 0.037596 | **0.010029** | **0.014399** | 0.115240 | 0.065660 | 0.182516 |

The strict holdout selected the ePSF (`1.428949`) over the independent
Gaussian (`1.445313`), radial profile (`1.469438`), and shared Gaussian
(`1.649514`). It recovered 25 sources after only 17 passed the circular seed
model. On the matched circular fixture, subpixel refinement also makes the ePSF
win narrowly (`1.421883` versus `1.423399`) and improve the final photometry and
shape rather than merely adding degrees of freedom.

RawTherapee 5.13 LMMSE and AMaZE were regenerated from this exact anisotropic
DNG with the checked-in neutral profiles. The ePSF wins aperture flux, FWHM,
elongation, RGB RMSE, and chroma MAE, while LMMSE wins stellar luminance and
false color. This establishes the intended construction-level advantage on one
controlled anisotropic fixture, not general state of the art. Multiple unseen
seeds, spatially varying PSFs, real RAWs, and another external implementation
remain necessary for a defensible general ranking.

## Low-Rank Chromatic ePSF Proof Of Concept

The chromatic extension retains one shared 2D ePSF and adds only a global X/Y
shift and isotropic scale for R and B relative to G. Scaling includes the area
Jacobian so source amplitudes remain flux parameters. Selection uses the same
unseen CFA samples plus a BIC penalty for the six additional parameters; this
rejects small apparent improvements on an achromatic low-information fixture.

The controlled fixture adds opposite `0.1875`-pixel R/B displacement and
opposite `0.10` relative scale around G to a mildly elliptical common Moffat PSF:

```bash
astrocfa-nogui benchmark-debayer --width 512 --height 384 --seed 7 \
  --noise astro --stars 25 --common-star-sigma 0.65 --moffat-beta 2.5 \
  --psf-ellipticity 0.20 --psf-angle 0.55 \
  --chromatic-psf-shift 0.1875 --chromatic-psf-scale 0.10 \
  --star-flux-scale 0.5 --star-phase diverse
```

The complete table, including RawTherapee outputs regenerated from this exact
DNG, is [`astrocfa-epsf-chromatic.csv`](results/astrocfa-epsf-chromatic.csv).

| Method | RGB RMSE | Chroma MAE | Star false color | Star luma RMSE | Flux relative | FWHM relative | Elongation |
|---|---:|---:|---:|---:|---:|---:|---:|
| Achromatic ePSF | 0.001478 | 0.000931 | 0.020008 | 0.018740 | 0.076537 | 0.010798 | **0.069004** |
| Chromatic ePSF | **0.001226** | **0.000922** | **0.013926** | **0.017387** | **0.072945** | **0.009229** | 0.071822 |
| RawTherapee LMMSE | 0.037607 | 0.002763 | 0.017528 | 0.019858 | 0.121535 | 0.071208 | 0.135041 |
| RawTherapee AMaZE | 0.038310 | 0.004518 | 0.018608 | 0.018655 | 0.125463 | 0.065731 | 0.148695 |

The penalized holdout selects the chromatic model (`1.423162 + 0.001618`)
over the achromatic ePSF (`1.427300`). It reduces false color by 30% relative to
its own ablation and beats both external methods in every reported stellar
metric; the achromatic ablation retains slightly better elongation. The recovered
R/B transform signs agree with the injected lateral displacement, although the
shared base ePSF absorbs part of the injected scale and Y shift.

This closes the observed color/luminance gap on one controlled fixture. It does
not establish general superiority: the chromatic model still needs multiple
unseen seeds, varying star colors and SNR, spatially varying chromatic aberration,
and real-camera validation.
