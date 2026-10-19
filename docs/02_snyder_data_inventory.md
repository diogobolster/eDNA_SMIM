# The Snyder flowing-water datasets (New Data folder)

Both are **steady-state additions**: eDNA release solution dripped at 100 mL/min with NaCl as a
conservative tracer, sampling after conductivity plateaued, n = 5 replicates per station,
sequential filtration through 1.2 um then 0.4 um polycarbonate, ddPCR (3 technical replicates
merged). Steady state and single fixed source are exactly the app's assumptions, and the design
matches Shogren et al. 2017. The two experiments are Chapters 5 and 6 of
`SNYDERE92024D.pdf` (Snyder 2024, PhD, Notre Dame).

## A. LEEF 2021, species-specific removal (Snyder Ch. 5)

File: `LEEF species specific eDNA removal eDNA concentrations and metadata.xlsx`

Design: four ND-LEEF streams (0.4 m wide, 0.06 m deep, 50 m, concrete-lined, ~1 L/s), each with
one substrate: Cobble (COB), Pea Gravel (PG), Sand (SAND), 3-way Mix (MIX). Four 3-week phases:
S = Early Shaded, LB = Early Open, HB = Late Open, RS = Late Shaded. At the end of each phase,
Common Carp (CARP, zoo pond water, 23 NTU) and Steelhead (STLHD, hatchery raceway, <1 NTU) eDNA
were added **simultaneously** to each stream. Stations A-E at 10, 20, 30, 40, 50 m.

So: 4 streams x 4 phases x 2 species = 32 additions, each with 5 stations x 5 reps x 2 sizes.

Sheets:
- `LEEF2021_BothSpeciesTotal` (776 rows): `GC.mL` = sum of the two size fractions per sample.
  Where a fraction was a non-detect it was summed as 0 (14 rows are 0 where both were missing).
- `LEEF2021_PltData_CARPSize_ForBa` (758) and `..._STLHDSize_ForB` (758): one row per sample and
  size. ddPCR raw fields present (`CopiesPer20uLWell`, `Positives`, `Negatives`, `AcceptedDroplets`,
  `Threshold` 6500 CARP / 3500 STLHD, `mLfiltered` 212-300). `GC.mL` is NA when copies = 0; the
  smallest quantified value has 2 positive droplets, and the largest non-detect also has 2 positive
  droplets, so the detection rule is on copies not droplets. Non-detects: **107 of 758 CARP rows,
  14 of 758 STLHD rows.** CARP 0.4 um concentrations have a median of 0.16 GC/mL, within a factor of 5 of the
  floor, and 73 of the 107 CARP non-detects are in the 0.4 um fraction.
- Hydraulics per stream and phase: `Q` (m^3/day), `w` (m) in the data sheets; and in the
  covariables file `U`, `D`, `Lambda`, `Beta` with SEs from rhodamine breakthrough curves, plus
  AFDM, chl-a, FBOM, GPP, ER. **These are the conservative SMIM parameters Eq. 4.15 needs**, for
  all 16 stream-phase combinations. `t1`, `t2` are not in the file (Volponi et al., in review).
- Size meaning: 1.2 um filter = particles > 1.2 um; 0.4 um filter = 0.4-1.2 um.
- Concentrations are raw GC/mL, not conductivity-corrected (concrete channel, no groundwater
  exchange, so no dilution along the reach).
- No source (bucket) rows: the `PltvBKD` column is "Plt" throughout. The thesis says release-
  solution samples were taken before and after each addition and reports their PSD by phase
  (Fig. 5.2E-F: CARP large > small in all phases but Early Open; STLHD sizes not different).
  **The bucket concentrations by size are the missing input for the Prediction tab.**

What Snyder found (Bayesian hierarchical fit, normalised to station 1, Vdep = kQ/w): total eDNA
median Sw 52 m CARP (16-515), 127 m STLHD (10-2174); Vdep 0.6-4.5 m/d CARP, 0.5-2.6 STLHD; no
clear species, biofilm or shade effect on total; large particles removed faster; the biofilm x size
interaction was significant for CARP only.

## B. Fort Custer, seasonal (Snyder Ch. 6)

File: `Fort Custer seasonal eDNA concentration data & metadata.xlsx`

Design: a natural first-order stream (Lawler Pond outflow, Fort Custer Training Center, MI; ~1 m
wide, 0.04 m deep, sand and gravel, forested). Steelhead only (absent from the watershed).
Four additions: F22 fall 2022, W23 winter, S23 spring, Sm23 summer 2023. Stations 20, 40, 60, 80,
100 m. Five reps per station. Fall was filtered 5/1.2/0.4 um; the 5 um fraction was folded into the
1.2 um fraction.

Sheets:
- `eDNApltdata_total` (100 rows), `eDNApltdata_bysize` (202 rows, 4 exact duplicate rows to drop).
- `plteDNA` = GC/mL **divided by the conductivity difference (plateau minus background)**. That is
  a dilution correction using the conservative tracer, done per sample. Units are therefore
  GC/mL per uS/cm, not GC/mL.
- Non-detects appear as NA: 8 of 25 in spring total (stations 60-100), 1 in winter at 100 m.
- Covariables file: per season Q (L/s), w, V (m/s), temperature, turbidity, DOC, pH, nutrients.
  Note the data sheet carries Q in m^3/day (91.584 = 1.06 L/s). No D, Lambda, beta here, so SMIM
  inversion is not possible for Fort Custer without a tracer breakthrough curve.
- No source rows. Thesis: source PSD by season in Fig. A.5 (summer large > small, p < 0.001; fall and
  winter marginal; spring no difference, indoor tank). Bucket samples were taken at start and end.

What Snyder found: total Sw 100 m summer, 125 fall, 179 spring, 227 winter; Vdep 0.53, 0.73, 0.72,
0.33 m/d (spring, summer, fall, winter); large ~2x faster than small in every season; composition
changed measurably with distance **only in winter** (beta regression, p < 0.001; p > 0.3 otherwise).
That last point is the single most important fact for the Prediction tab: if composition does not
change along the reach, the method has nothing to work with.

## C. Also in the folder but not flowing-water

`doi_10_5061_dryad_37pvmcvnw__v20221219/`: the polydisperse mesocosm time series behind Section 4.2
(Brandão-Dias et al. 2023a). `es3c02638.pdf` is the ES&T suspended-materials paper whose data is
the app's demo template.

## D. Step-0 baseline, app regression applied as-is (`03_baseline_fits.csv`)

Pooled-replicate semi-log regression per dataset and component, non-detects dropped, lambda
clipped at 0, exactly `fit_one()`. Headline:

| | fits | lambda 95% CI excludes 0 |
|---|---|---|
| LEEF 1.2 um | 32 | 13 |
| LEEF 0.4 um | 32 | 13 |
| LEEF total | 32 | 11 |
| Fort Custer 1.2 um | 4 | 4 |
| Fort Custer 0.4 um | 4 | 2 |
| Fort Custer total | 4 | 3 |

- Fort Custer reproduces Snyder's totals well without any Bayesian machinery: Vdep 0.53 / 0.80 /
  0.60 / 0.24 m/d against her 0.53 / 0.73 / 0.72 / 0.33; Sw 170 / 100 / 154 / 194 m against
  179 / 100 / 125 / 227. By size, 1.2 um lambda > 0.4 um lambda in all four seasons, as she found.
- LEEF is noisy: R^2 is mostly < 0.4, several profiles are flat or rise with distance (Sand and Mix
  streams especially; e.g. CARP S_SAND 1.2 um has a significantly *positive* slope), and many 0.4 um
  CARP series lose a third of their points to non-detects. Where the fit is significant, Sw is
  17-88 m, in line with her medians (52-71 m by size).
- The two sizes have different lambda in the same addition often enough to give the composition a
  trend, but the differences are small (typically 0.01-0.02 m^-1) relative to the CI widths.
