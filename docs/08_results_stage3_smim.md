# Results: Stage 3, SMIM-rates tab on LEEF 2021 (as-is method)

Code: `edna_py/run_stage3_smim.py`. Outputs: `results/stage3_rs_from_mesocosm.csv`,
`results/stage3_smim_rh.csv`. Eq. 4.15 exactly as in `appv6+.R` (`k_func`), t1 = 1 s,
t2 = 1e5 s as hard-coded there, r_h found by least squares on |k_model - k_data|.

## Inputs

k_data: minus the slope of the app's regression (Stage 2) per stream x phase x species x
{1.2 um, 0.4 um, total}: 96 combinations. 75 have k > 0; 40 have a slope CI excluding 0.
Only those 40 are reported as results; the rest are in the CSV, flagged.

Conservative parameters: U, D, Lambda, beta per stream and phase from Snyder's rhodamine fits
(`LEEF2021Ancillary`): U 0.06-0.20 m/s, D 0.010-0.039 m^2/s, Lambda 0.015-0.085 1/s,
beta 0.64-1.10.

r_s, the water-column removal rate, fitted here from the Dryad polydisperse mesocosm series
(Brandao-Dias 2023a, treatment AF, 1/6/18/30 h, semi-log regression pooling replicates):

| species | component | Dryad fractions used | r_s (1/h) | 95% CI | R^2 | r_s (1/s) |
|---|---|---|---|---|---|---|
| CARP | 1.2 | 10 + 1 um summed | 0.081 | 0.060-0.102 | 0.88 | 2.2e-5 |
| CARP | 0.4 | 0.4 um | 0.053 | 0.032-0.075 | 0.47 | 1.5e-5 |
| CARP | total | unfiltered (N) | 0.090 | 0.075-0.104 | 0.82 | 2.5e-5 |
| STLHD | 1.2 | 10 + 1 um summed | 0.077 | 0.034-0.120 | 0.61 | 2.1e-5 |
| STLHD | 0.4 | 0.4 um | 0.131 | 0.109-0.154 | 0.84 | 3.7e-5 |
| STLHD | total | unfiltered (N) | 0.091 | 0.077-0.104 | 0.85 | 2.5e-5 |

These sit at the app's default r_s of 3e-5 1/s. The size mapping (Dryad 10 + 1 um filters for
LEEF's >1.2 um, Dryad 0.4 um for LEEF's 0.4-1.2 um) is approximate; it does not matter, see the
sensitivity below.

## Results

For the 40 significant profiles, r_h in 1/s:

| species | component | n | median r_h | range |
|---|---|---|---|---|
| CARP | 1.2 | 7 | 6.8e-3 | 4.2e-3 to 1.4e-2 |
| CARP | 0.4 | 4 | 7.3e-3 | 3.3e-3 to 2.9e-2 |
| CARP | total | 7 | 4.3e-3 | 2.7e-3 to 7.5e-3 |
| STLHD | 1.2 | 7 | 4.7e-3 | 2.3e-3 to 1.1e-2 |
| STLHD | 0.4 | 8 | 3.0e-3 | 1.7e-3 to 1.5e-2 |
| STLHD | total | 7 | 4.0e-3 | 1.3e-3 to 7.1e-3 |

Per-profile ranges from the k CI are a factor 3-10 wide (column `rh_at_k_CI_low/high`).

Three points.

1. **r_h / r_s is 50-600, median 190.** The hyporheic rate is two orders of magnitude above
   the water-column rate, the same ratio Sec. 4.4 found for ARGs. Mechanically this is forced:
   the observed in-stream loss rate k U is about 2e-3 1/s (median) against a water-column
   decay r_s of 2.5e-5 1/s, a factor 80. With r_h = 0 the model predicts k of 1e-4 to 6e-4
   1/m, 50-100 times below the measured 0.01-0.05 1/m. The SMIM has exactly two sinks, water
   column and hyporheic zone, so whatever removes eDNA that is not water-column decay is
   booked as r_h. Deposition onto the bed, biofilm capture and adsorption to benthic organic
   matter, which Snyder's Vdep framing attributes the removal to, all land in r_h here.
   The number is therefore "removal rate while in storage", not a biological decay rate.
2. **r_s is irrelevant and t2 matters a little.** Setting r_s = 0 changes r_h by 1.6%
   (median); doubling it, even less. t2 = 1e3 s instead of 1e5 s raises r_h by 27% (median).
   The rebuild's guess of 1000 s would have shifted every number by that amount; Daniel's
   1e5 s is what the app ships, and the true t2 from Volponi would settle it.
3. **r_h follows k, not the hydraulics.** Spearman r_h vs k_data = 0.70; vs Lambda -0.21;
   vs beta 0.06. The sand stream has the highest median r_h (1.4e-2 vs 4-5e-3 for the other
   three), driven by two steep profiles (CARP S_SAND 0.4 um, CARP HB_SAND 1.2 um) while most
   sand profiles are flat or rising and drop out. By phase: HB 3.0e-3, LB 4.1e-3, S 4.8e-3,
   RS 5.8e-3, no monotone biofilm trend. Size: CARP large > small is not supported (medians
   6.8e-3 vs 7.3e-3, n = 4 for small); STLHD large > small (4.7e-3 vs 3.0e-3), in line with
   Snyder's Vdep ordering.

21 of 96 profiles have k <= 0 (concentration flat or rising downstream); the SMIM cannot
represent those at all, since both sinks only remove. They are the sand and mix streams in
most phases, for both species.

## What this means for the write-up

The SMIM tab runs on LEEF and returns consistent numbers, but what it returns is a
hyporheic-zone removal rate that absorbs every non-water-column loss. Its headline, r_h two
orders of magnitude above r_s, is a restatement of "in-stream removal is 50-100x faster than
mesocosm decay", which Snyder's Vdep already says. The tab adds a process attribution (storage
zone) that the data cannot test, because deposition is not a separate term in Eq. 4.15.
Record it as such; do not present r_h as a measured hyporheic decay rate.

Without a t2 from the rhodamine fits the absolute r_h carries a ~30% convention uncertainty on
top of the factor 3-10 from the regression CI.
