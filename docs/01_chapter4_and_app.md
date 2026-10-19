# Chapter 4 (Hallack 2025) and the R Shiny app, as written

Source: `METANIAS CARVALHO HALLACKD42025D.pdf`, Chapter 4, pp. 68-110, and `edna-app/` (R) with
`edna-app-py/` (Python port). Nothing below goes beyond what the thesis and the app do.

Provenance. The `edna-app/` folder is a **reconstruction** from the thesis text and the live
interface, not Daniel's source. Daniel's real source is github.com/dhallackla/edna_app, now copied
to `Starting Anew/edna_app_original/`; `06_original_app_vs_rebuild.md` lists where the real code
departs from the thesis text and from the rebuild (replicate-max normalisation before regression,
n-1 fractions in the distance, boundary solutions discarded, t2 = 1e5 s, SANN optimiser,
histogram mode). The sections below describe the method as the thesis states it; the code note
takes precedence where they differ.

## 4.2 Which decay model (Brandão-Dias et al. 2023a, Mol. Ecol. Res.)

Four candidate models for n size classes: (1) one shared removal rate, (2) independent first-order
removal per class, (3) cascade (large breaks into small, only smallest removed), (4) cascade plus
independent removal. Fitted to the polydisperse mesocosm experiment (`New Data/doi_10_5061_dryad_37pvmcvnw`:
CARP and STHD, pre-filtered 10/1/0.4/0.2 um, sampled 1, 6, 18, 30 h, sequentially filtered). By AIC/BIC
the independent model (2) wins overall; model (1) is chosen least often. That is the justification for
everything that follows: **each component decays independently at its own first-order rate.**

## 4.3 The source-location method (the part the app implements)

Setting. A fixed source releasing continuously at rate p into discharge Q, so C_source = p/Q. Steady
state. Single source. Reference x = 0 is the most upstream sample; the source is at x < 0.

Eq. 4.5    C_i(x) = C0_i exp(-lambda_i x),  i = 1..n;   C(x) = sum_i C_i(x)

Eq. 4.6    ecd(x) = (f_1, ..., f_{n-1}),  f_i = C_i / C   (the component distribution; n-1 free dims)

Eq. 4.8    accept a trajectory if  min_x ||ecd(x) - ecd_source|| < eps,  x* = argmin

Estimation (Algorithm 2):
1. For each component, semi-log regression log C_i vs x on the sampling data. Take lambda_i = -slope,
   C0_i = exp(intercept), and their 95% confidence intervals. lambda is truncated at 0 (a single
   source cannot make concentrations rise downstream).
2. Draw N_trials parameter sets, each lambda*_i ~ U(lambda_low, lambda_up), C0*_i ~ U(C0_low, C0_up),
   independently per component.
3. Back-extrapolate each set over x in [-x_lim, 0] (C_i(-d) = C0_i exp(+lambda_i d)).
4. Form ecd(x) along each trajectory.
5. Keep trajectories that pass within eps (Euclidean distance in composition space) of ecd_source;
   record x* and C(x*).
6. Histograms of x* and C(x*) are the answer. Summaries: P10, P50, P90 and the mode of a lognormal
   with the same mu, sigma as log(x*).

Required inputs: component concentrations at >= 3 positions (or times); the source composition
ecd_source; its tolerance eps; x_lim, step, N_trials.

Validation in the thesis (4.3.1.2, 4.3.2). Recirculating mesocosms (Snyder et al. 2023 /
Brandão-Dias et al. 2025): 6 treatments x 4 reps, STHD, sampled 40 min, 6, 18, 48 h, filtered
10/1/0.2 um, two amplicon lengths. **The 40 min sample was declared the source** (ecd_source = mean
composition over 4 reps, eps = SD over reps), the other three times were the sampling data, shifted
so t = 6 h became 0. Known answer: 5.33 h. Settings: x_lim 100 h, step 0.5 h, 1000-2000 trials.
Result (Table 4.2, 3 components): [P10, P90] brackets 5.3 h in 9 of 12 datasets; in 10 of 12 the
mode is within 5 h of the truth; average [P10, P90] = [2, 44] h against a prior of [0, 100]; the
two Open Bare Substrate failures returned only 3-4 accepted trials, Open Biofilm Long was too
narrow ([0.5, 2] h), Shaded No Substrate Long too wide ([2, 88] h). More components gives narrower but more biased answers
(Table 4.3, Fig. 4.10). Stated limitations: 3 time points is thin; single fixed steady source;
independent decay; uniform independent parameter draws.

Two lines in the thesis that matter for the field data:
- "The nearest or most dominant source typically drives the signal" (p. 100): multiple sources are
  acknowledged, not handled.
- "Targeting different mitochondrial genes, which maintain a consistent 1:1 ratio ... provides an
  internal reference that eliminates uncertainty in initial proportions" (p. 98).

## 4.4 Reactive SMIM at steady state (ARG work, Badilla-Aguilar et al. 2024)

Eq. 4.15    C(x) = C0 exp(k x),   k = (V - sqrt(V^2 + 4 D Theta)) / (2D)
            Theta = r_s + Lambda [ 1 - (1 + t2 r_h)^beta e^{t1 r_h} Gamma(-beta, t1(1/t2 + r_h)) / Gamma(-beta, t1/t2) ]

Conservative parameters V, D, Lambda, beta, t1, t2 come from a rhodamine breakthrough curve. k comes
from the semi-log regression of the concentration profile. r_s (water-column removal) is taken from
mesocosm experiments. r_h (hyporheic removal) is then the only unknown, found by least squares.
Finding in the ARG data: r_h is two orders of magnitude above r_s.

## 4.5 The app

Inputs: a sampling CSV (`Treatment, Replicate, Filter, Target, Concentration, Distance|Time`) and a
source CSV (`Treatment, Replicate, Filter, Target, Concentration`). Data type is "total" or "by size";
the x axis is distance or time. Regressions run automatically on upload.

| Tab | Needs | Gives |
|---|---|---|
| LOD-distance | total or by-size vs x | x_LOD = ln(C0/LOD)/lambda, with a range from the CI corners |
| SMIM-rates | by-size vs **distance** + V, D, Lambda, beta, t1, t2, and one of r_s/r_h | the other rate |
| Prediction | by-size + source composition | histograms of x* and C(x*) (Algorithm 2) |
| Sensitivity | by-size | composition upstream for chosen x and rates (forward only) |

Defaults recovered from the live app: nsim 1000, tmax 100, step 2, tol 0.05, LOD 0.05 copies/mL,
U 0.0982 m/s, D 0.0345 m^2/s, Lambda 0.085, beta 0.898, r_s 3e-5, r_h 1e-3.

Stated assumptions (p. 110): independent first-order decay per category; fixed source location;
steady state; conservative SMIM parameters known.

The demo sampling template is the ES&T 2023 suspended-materials experiment (treatments A_NP, B_PG,
C_FPOM, D_Algae, E_Clay, F_TiO2 x CARP/STHD; "Distance" = 0, 5, 17, 29 are hours after the first
sample). The source template is the bucket composition for the same experiment.

## How the regression is actually done in the code (matters for step 0)

`fit_one()` in `edna-app/R/model.R`: `lm(log(conc) ~ x)` on **all replicate points pooled**, dropping
conc <= 0 and NA; `confint()` at 95%; lambda = max(-slope, 0), bounds likewise clipped at 0;
C0 = exp(intercept), bounds from the intercept CI. Needs >= 3 points and >= 2 distinct x.
Non-detects are therefore simply discarded. The LOD tab regresses the per-replicate total
(components summed within a replicate at each position).
