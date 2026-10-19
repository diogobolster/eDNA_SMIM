# Daniel's original app (github.com/dhallackla/edna_app) vs the rebuild in `../edna-app`

Copied without git history to `edna_app_original/` (commit 4ac13a7, 21 May 2024). Main file
`appv6+.R` (3200 lines), model helpers in `aux_fun.R`, help text in `help_mds/`, five CSV
templates. Packages: shiny, shinyWidgets, tidyverse, DescTools, optimx, gt, broom, shinyhelper,
ggforce, markdown, ggnewscale, expint. **This replaces the rebuild as the reference. Use it.**

What the original actually does, where it differs from the thesis text or the rebuild:

## Regression (both tabs)

- Concentrations are first **normalised within each Treatment x Replicate by that replicate's
  maximum** (over all filters and all distances). The regression is `lm(ln_norm_c ~ Distance)`
  per Filter, pooling replicates. Zeros become NA and are dropped. CIs are `broom::tidy(conf.int
  = TRUE)`, 95%. Only rows with Distance >= 0 are used.
- Consequence: slopes are those of the normalised series; because every point in a replicate is
  divided by one number, the slope within a replicate is unchanged, but the pooled fit differs
  from a pooled fit on raw data whenever replicates differ in level. Composition fractions are
  unaffected (the normaliser cancels within a sample).
- The slope is **not** truncated in the regression table. Truncation happens in the Monte Carlo:
  each drawn slope is clipped with `min(slope, 0)` (i.e. lambda >= 0) in `single_pred()`.
- Denormalisation back to copies/mL uses the mean of the replicate maxima (`denorm_factor`).
- The rebuild regresses raw log concentrations and clips lambda and its CI at 0 in the table.
  `03_baseline_fits.csv` was computed the rebuild's way; Stage 2 must redo it the original's way.

## Prediction (Algorithm 2)

- `MC_new()`: `set.seed(1)`; for each trial draws intercept and slope for every filter
  independently and uniformly within their 95% CIs; back-extrapolates over `time = -seq(0, tmax,
  dt)`; composition from the normalised concentrations.
- Distance to the source composition is Euclidean over the **first n-1 fractions**
  (`x[3:(Nsizes+1)]` vs `f_initial$mean[1:(Nsizes-1)]`), matching Eq. 4.6's n-1 dimensions.
  The rebuild uses all n fractions, which scales the distance by sqrt(2) for n = 2.
- Acceptance: per trial, the grid point of minimum distance; accept if `<= tol`; then
  **solutions at the two boundaries, x = 0 and x = -tmax, are discarded**
  (`filter(time != -tmax, time != 0)`). This is not in the thesis or the rebuild and it matters:
  a trajectory whose closest approach is at the edge of the search window is treated as no
  solution. Optional cutoffs `d_max` and `c0_max` as in the rebuild.
- Source composition from the source CSV: fraction per Treatment x Replicate, then mean and SD per
  Filter across replicates. The SD is displayed; `tol` is a manual input, default 0.05.
- Summaries shown by the app: mean, SD and a **histogram mode (15 equal bins)** for distance and
  for log10(C0). The P10/P50/P90 of Tables 4.2-4.3 are not produced by the app; they were
  computed from the accepted trials outside it. The thesis "mode" is this histogram mode, not
  exp(mu - sigma^2) as the rebuild assumed. Lognormal or normal curves are overlaid on the
  histograms for display only.
- Defaults: nsim 1000, tmax 100, dt 2, tol 0.05 (as the rebuild found from the live app).

## LOD-distance

`x_LOD = (ln LOD - intercept - mean(ln max_c)) / slope`, with an upper value from the CI corner
(low intercept, high slope) and `max(ln max_c)`, and a lower value from the opposite corner and
`min(ln max_c)`. Negative results are NaN. Per filter when data are by size. The rebuild's
formula is algebraically the same once the normalisation is undone, but its range uses only the
regression CI, not the spread of replicate maxima.

## SMIM-rates

- `k_func()` is Eq. 4.15 with `expint::gammainc`; k is returned positive and used as
  `C = C0 exp(-k x)`.
- **t1 = 1 s and t2 = 1e5 s are hard-coded** (`logT1 <- 0; logT2 <- 5`). The rebuild guessed
  t2 = 1000 s. The UI inputs for them are commented out.
- The unknown r_h is found by `optim(method = "SANN")` (simulated annealing) minimising
  `|k_model - k_data|`, with k_data = -slope of the (normalised) regression for the chosen filter,
  and r_s from the input (default 3e-5). Starting value 1e-3. SANN is stochastic and has no
  convergence check; the rebuild's 1-D `optimize` on log r_h is the better numerical method for
  the same objective, but the reference answer is whatever SANN returns.
- Defaults U 0.0982, D 0.0345, Lambda 0.085, beta 0.898.

## Templates

- `template_sampling_perSize.csv` and `template_source.csv` are byte-identical to the rebuild's
  demo files (the ES&T suspended-materials mesocosm experiment).
- `template_sampling_Total.csv` is **LEEF 2021 total eDNA**: 16 treatments `HB_<sub>_<sp>` and
  `SH_<sub>_<sp>`, reps R1-R4, Distance 0-40 (station minus 10 m). Daniel had this data in May
  2024, so the LOD tab has already been pointed at it. The numbers are **not** the ones in
  Snyder's final file (different values, 27 NAs, 4 reps not 5, two phases not four): a
  preliminary cut. Use Snyder's file, not this template.
- `*_hours.csv` variants carry an `Hours` column instead of `Distance` for time-series data.

## What this changes in the plan

- Decision 3 is settled: the original is the engine. Stage 0 runs `appv6+.R` (needs R with the
  twelve packages above; the rebuild README's conda recipe works, adding the extra packages).
- Stage 2 recomputes the baseline with the original's max-normalisation and boundary-exclusion
  rules so the numbers match the app exactly.
- Stage 3 uses t1 = 1 s, t2 = 1e5 s as the app's reference values and still asks Volponi for
  the measured ones.
- The boundary-exclusion rule in Prediction should be recorded as part of the as-is method; it
  will remove some of the n = 2 trajectories that would otherwise be accepted at x = 0.
