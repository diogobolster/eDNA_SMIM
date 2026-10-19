# Results: Stages 0, 1, 2 and 4a (as-is method, Python port of Daniel's code)

Code: `edna_py/original.py` (the port), `edna_py/run_stage0_template.py`, `edna_py/run_snyder.py`.
Outputs: `results/`. Everything below is the published method run unchanged; nothing has been
altered or tuned. Distances in the tables are reported from the injection point unless the
column says otherwise.

## Stage 0. The port and the template

Port checks against R. The regression (lm + 95% t-based CIs), the upper incomplete gamma for
negative order, Eq. 4.15's k, and the r_h inversion all reproduce the R reference values in
`../edna-app-py/tests/r_reference.tsv` to 12 significant figures. The Monte Carlo and
solution-selection logic is a line-by-line transcription of `MC_new()`, `single_pred()`,
`df_close()` and `df_close_C0()`; it cannot be checked bit-for-bit because R's `set.seed(1)`
stream is not reproducible in numpy. SANN is replaced by a bounded deterministic search on the
same objective.

Template (ES&T 2023 mesocosms, 12 treatments, 3 size classes, known answer 1 h before the
first sample), thesis settings 2000 trials, 0.5 h step, tol 0.05:
`results/stage0_template_prediction.csv`.

| | treatments |
|---|---|
| no regression possible for 2 of 3 classes (clay, TiO2: eDNA all on the 10 um filter) | 4 |
| fewer than 4 accepted trials | 3 |
| [P10, P90] contains 1 h at tol 0.05 | 0 of 12 |
| [P10, P90] contains 1 h at tol = source SD | 2 of 12 |

Where there are solutions the P10 sits at 1.5-4 h and the mode at 2-5 h, so the answer is
missed by a factor 2-4 on a 0.5 h grid. This demo is a harder case than the thesis
validation: the bucket composition is not the water-column composition once clay, TiO2 or FPOM
adsorb the DNA (the paper's own finding), and 1 h is one grid step from the excluded x = 0
boundary. It is a check that the pipeline runs end to end, not a validation; the thesis
validation data (Table 4.2) are not in the folder.

## Stage 1. App-format inputs

`results/app_inputs/{LEEF2021,FortCuster}_sampling_{perSize,Total}.csv`. Distance = 0 at the
first station (app convention). Non-detects written as 0 (the app's convention; it drops them
from the regression and keeps them in the composition). These load into `appv6+.R` unchanged.

## Stage 2. Regressions and LOD, the original's way

`results/stage2_regressions_original_method.csv` (108 fits: 36 datasets x {1.2, 0.4, total}).

The replicate-max normalisation changes lambda by a median of 0.0007 m^-1 against the raw fits
in `03_baseline_fits.csv` (max 0.066 in one flat LEEF profile) and moves three fits across the
significance line (46 -> 49 of 108 with the slope CI excluding 0). The two conventions agree
in substance.

Fort Custer, original method (lambda in m^-1, Sw in m, Vdep in m/d; CI = slope 95% CI
excludes 0):

| season | component | lambda | R^2 | Sw | Vdep | CI excl. 0 |
|---|---|---|---|---|---|---|
| F22 | 1.2 | 0.0094 | 0.17 | 106 | 0.87 | yes |
| F22 | 0.4 | 0.0034 | 0.04 | 296 | 0.31 | no |
| F22 | total | 0.0065 | 0.20 | 154 | 0.60 | yes |
| S23 | 1.2 | 0.0095 | 0.21 | 105 | 0.86 | yes |
| S23 | 0.4 | -0.0041 | 0.04 | inf | | no (rising) |
| S23 | total | 0.0057 | 0.12 | 176 | 0.51 | no |
| Sm23 | 1.2 | 0.0109 | 0.28 | 92 | 0.87 | yes |
| Sm23 | 0.4 | 0.0069 | 0.38 | 144 | 0.56 | yes |
| Sm23 | total | 0.0100 | 0.40 | 100 | 0.80 | yes |
| W23 | 1.2 | 0.0122 | 0.59 | 82 | 0.57 | yes |
| W23 | 0.4 | 0.0068 | 0.31 | 146 | 0.32 | yes |
| W23 | total | 0.0054 | 0.26 | 185 | 0.25 | yes |

Against Snyder (Bayesian, normalised to station 1): total Vdep 0.51/0.80/0.60/0.25 here vs
0.53/0.73/0.72/0.33 (spring/summer/fall/winter); 1.2 um Sw 82-106 m here vs her 55-69 m; the
large fraction decays faster than the small in every season, as she found. The app's plain
regression gets the same story.

LEEF: 14 of 32 total fits, 14 of 32 at 1.2 um and 12 of 32 at 0.4 um have a slope CI excluding
0. Where significant, Sw is 17-90 m. The Sand and Mix streams are mostly flat.

LOD distances (`results/stage2_lod_distances.csv`). For LEEF at the ddPCR floor (0.035 GC/mL
CARP, 0.056 STLHD) the predicted detection limit on total eDNA is a median 205 m (CARP,
range 81-945) and 248 m (STLHD, 153-1483) from the injection point, i.e. well beyond the 50 m
reach; the upper/lower bounds from the CI corners are often undefined because one corner has
a non-negative slope. **For Fort Custer the LOD tab is not interpretable**: `plteDNA` is GC/mL
divided by the conductivity difference, so its values (0.005-0.06) are not in copies/mL and no
LOD in copies/mL applies. Undoing that needs the conductivities, which we will not get.

## Stage 4a. Reference-station prediction

`results/stage4a_reference_station_prediction.csv`. Source = composition at the first station
(mean over replicates with both fractions quantified; SD recorded), sampling data = the other
four stations re-zeroed at the second, app defaults (1000 trials, 100 m window, 2 m step,
tol 0.05). Known answer: 10 m (LEEF), 20 m (Fort Custer). 34 of 36 datasets could be run;
CARP LB_COB (too few 1.2 um detections once station 1 is removed) and CARP RS_PG (no station-1
replicate with both fractions quantified) could not.

| | datasets | coverage [P10,P90] | median P10 | median P50 | median P90 | median rel. width | median n accepted |
|---|---|---|---|---|---|---|---|
| LEEF 2021 | 30 | 24 of 30 (0.80) | 5 | 20 | 66 | 0.55 | 316 |
| Fort Custer | 4 | 3 of 4 (0.75) | 14 | 49 | 82 | 0.68 | 264 |

Fort Custer in full (x = metres upstream of the second station; truth 20):

| season | f_0.4 at source | n accepted | P10 | P50 | P90 | mode | covered |
|---|---|---|---|---|---|---|---|
| F22 | 0.40 | 247 | 10 | 48 | 82 | 37 | yes |
| S23 | 0.18 | 99 | 56 | 80 | 94 | 88 | no |
| Sm23 | 0.26 | 408 | 6 | 28 | 78 | 5 | yes |
| W23 | 0.48 | 281 | 18 | 50 | 82 | 37 | yes |

Three things the as-is run shows, all properties of the method rather than of this port:

1. **With two components the tolerance does nothing.** Results at tol 0.05 and at tol = source
   SD are identical in every dataset. With n = 2 the composition is one number, a trajectory
   is monotonic in it, so its closest approach to f_source is either the grid point nearest
   the crossing (residual at most half the composition change over one grid step, far below
   either tolerance) or an endpoint of the window (discarded by the boundary rule). The
   tolerance only matters if it is set below that grid resolution, and then it removes
   trajectories whose crossing falls between grid points without moving any solution
   (`tests/test_original.py::test_two_components_make_tolerance_irrelevant`). The prediction
   reduces to "where do the sampled (intercept, slope) pairs make the trajectory cross
   f_source", and the spread comes only from the regression CIs.
2. **The histograms are front-loaded.** In 22 of 34 datasets the histogram mode is the first
   bin (x of 3-6 m), so the relative error of the mode is -0.48 almost everywhere and P50 is
   several times the mode. This is the asymmetry the thesis reports (mode << P50, Sec. 4.3.2);
   here it is extreme, and the mode is not a useful point estimate.
3. **Coverage is bought with width.** [P10, P90] spans 55-68% of the 100 m prior window
   (thesis case D, also two components: 0.58). The method knows the source is upstream and
   not very far, and little else. Fort Custer spring, the one season with a rising 0.4 um
   profile, puts the source 56-94 m upstream when it is 20 m.

Whether the composition changes with distance at all (Snyder's beta regression: only in
winter) does not separate the covered from the uncovered datasets here; the LEEF datasets that
miss have slightly larger lambda contrasts (median 0.014 vs 0.010 m^-1) but fewer accepted
trials. The misses are the datasets where the first-station composition is far from where the
back-extrapolated trajectories go, i.e. where the between-station composition is not
monotone.

## What this leaves for Stage 3 and the write-up

- SMIM inversion on LEEF needs r_s per size. The Dryad polydisperse time series in the folder
  gives water-column decay rates by filter size for both species; those are the natural r_s.
  t1 = 1 s, t2 = 1e5 s as in the app.
- No bucket data means Stage 4b cannot run; 4a is the field test of the method and is done.
- The write-up should say plainly: on 36 field datasets with two size classes, the as-is method
  returns an 80% interval that covers the truth four times in five and is roughly half the prior
  window wide; the mode is uninformative; tolerance is irrelevant at n = 2.
