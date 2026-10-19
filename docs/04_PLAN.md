# Plan: apply Chapter 4 and the app, unchanged, to the Snyder flowing-water data

Rule for this folder: the model is the thesis model and the app is the app. Every stage runs the
published method as-is and records what it gives. Changes to the method are a separate, later
decision, made only after the as-is results are on the table.

## Stage 0. Prove the tool on its own demo (half a day)

Engine: **Daniel's original `appv6+.R`** (github.com/dhallackla/edna_app, copied to
`edna_app_original/`; see `06_original_app_vs_rebuild.md` for how it differs from the thesis
text and from the rebuild). Run it on its bundled template. The template is the ES&T
suspended-materials experiment with the bucket composition as source and "Distance" = hours after
the 1 h sample, so the known answer for the Prediction tab is 1 h. Check that the regressions,
LOD and Prediction outputs are sensible and that the thesis's settings (1000-2000 trials, tol
from the replicate SD) reproduce a histogram that brackets 1 h. This is the only validation with a
known answer we can run without more data. If the Snyder et al. 2023 mesocosm data (the 24
recirculating mesocosms behind Table 4.2) can be obtained, reproduce Table 4.2 for a few rows too.

## Stage 1. Put the Snyder data into the app's format (half a day)

One script, `make_app_inputs.py`, that writes app-format sampling CSVs:
`Treatment, Replicate, Filter, Target, Concentration, Distance`.
- LEEF: Treatment = `<Phase>_<Substrate>` (16), Target = CARP or STLHD, Filter = 1.2 or 0.4,
  Replicate = R1-R5, Distance = station minus 10 m so the first station is x = 0 as the thesis
  prescribes. 32 datasets, two components each. Non-detects written as 0 (the app drops them).
- Fort Custer: Treatment = season (4), Target = STHD, Filter = 1.2 or 0.4, Distance = station
  minus 20 m. Concentration = `plteDNA` (conductivity-corrected). 4 datasets.
Keep a second column set with the raw station distances so results can be reported from the
injection point.

## Stage 2. LOD-distance tab on all 36 datasets (one day)

This needs nothing but the data. For every dataset and for total, 1.2 um and 0.4 um: lambda, C0,
95% CIs, R^2, x_LOD at the app's default LOD and at the ddPCR floor actually seen in the data
(about 0.035 GC/mL for CARP, 0.056 for STLHD at 250 mL filtered). `03_baseline_fits.csv` is a first pass done
the rebuild's way (raw concentrations); the original app regresses replicate-max-normalised
concentrations, so the stage is to recompute with that rule, run the same data through the app,
confirm the numbers match, and tabulate against Snyder's Sw and Vdep. Expected outcome: Fort Custer agrees with her;
LEEF agrees where her fits were good and disagrees where the profile is flat, because her
normalisation to station 1 and her priors pull k positive and the app's truncation at 0 does the
same thing more crudely.

## Stage 3. SMIM-rates tab on LEEF (one to two days)

LEEF is the one dataset with the conservative SMIM parameters measured: U, D, Lambda, beta per
stream and phase (16 sets). For each stream-phase-species-size (64 combinations, fewer after
dropping flat profiles): k from Stage 2, r_s from the mesocosm removal rates by size (Snyder
Ch. 2/3, Brandão-Dias 2023a), invert Eq. 4.15 for r_h. Two inputs are missing and have to be
obtained before this stage is trustworthy: t1 and t2 (the original app hard-codes 1 s and 1e5 s;
Volponi et al. have the measured values) and the mesocosm r_s by size in s^-1. Output: a table of r_h by
species, size, substrate and biofilm phase, with the sensitivity of r_h to t1, t2 and r_s.
Nobody has computed hyporheic removal rates for eDNA by particle size; this is a result in its own
right if the inversion is stable.

## Stage 4. Prediction tab (two to three days)

The thesis validated the method by declaring the first time point the source and predicting its
distance from the later points (Sec. 4.3.1.2). Do the identical thing in space:

4a. **Reference-station test, no new data needed.** Source = mean composition at the first
station (10 m LEEF, 20 m Fort Custer), eps = SD of composition across its 5 reps (as in Table 4.1).
Sampling data = the remaining four stations, shifted so the second station is x = 0. Known answer:
10 m (LEEF) or 20 m (Fort Custer). Thesis settings: x_lim 100 m, step 1-2 m, 1000-2000 trials.
Run all 36 datasets. Report P10/P50/P90, mode, number of accepted trials, and the three thesis
metrics (truth coverage, relative error of the mode, relative width).

4b. **True-source test.** Needs the release-solution (bucket) concentrations by size for each
addition from Elise (she reports they were taken before and after each addition; the by-phase
PSDs are in her Fig. 5.2E-F and Fig. A.5). Then source = bucket composition, data = all five
stations, known answer = 10 m or 20 m from the injection point. Also ask for t1, t2 and the
Fort Custer plateau/background conductivities so the correction can be undone if needed.

Expected outcome, stated now so it can be checked: the method only has information where the two
sizes decay at different rates. From Stage 2, that holds clearly in Fort Custer every season
(1.2 um lambda 1.6-2.8x the 0.4 um lambda) and rarely in LEEF: both size fits are significant in
only 4 of 32 additions, at least one in 22, and the median difference between the two lambdas
is 0.01 m^-1, smaller than most of the CI half-widths. Snyder's
own beta regression found a composition trend only in Fort Custer winter. So expect: a usable
answer for Fort Custer winter, maybe summer; wide or empty histograms for most of LEEF, like the
thesis's Open Bare Substrate case. With two components the composition is one number (f_1.2), so
the tolerance region is an interval on a line and every trajectory that moves at all will cross
it somewhere; the accepted-trial count will be high and the distance histogram wide. That is a
property of n = 2, not a bug, and it is the first thing the as-is run will show.

## Stage 5. Write up what the as-is method does on field data (one day)

A short results note: per-dataset table, the three thesis metrics, and a plain list of where it
worked, where it did not, and the data feature responsible in each failure (flat profile,
non-detects, two components, tolerance from replicate SD being wider than the whole trajectory).
Only at this point decide what, if anything, to change in the method. The candidates will be
obvious from the failures; do not pre-empt them.

## Decisions needed from Diogo before Stage 1

1. Engine is settled: the original `appv6+.R`. Open question is only whether to drive it as the
   Shiny UI (one dataset at a time, 36 datasets) or to source `aux_fun.R` and the regression
   block into a headless R script that loops over datasets with identical settings. The headless
   route is the practical one; it reuses Daniel's functions verbatim.
2. Send the data request to Elise now (bucket concentrations by size per addition; t1, t2; mesocosm
   r_s by size; Snyder 2023 mesocosm data for Table 4.2). Draft in `05_data_request.md`.
3. Distances measured from the first station (thesis convention; the app itself filters
   Distance >= 0) or from the injection point for reporting. The plan uses the thesis convention
   for fitting and reports both.
