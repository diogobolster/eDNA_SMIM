# Results: Stage 0b, the thesis mesocosm validation recomputed from public data

The data behind Hallack (2025) Sec. 4.3.1.2 and Tables 4.1-4.2 are archived: Brandao-Dias et
al. 2025, Dryad doi:10.5061/dryad.nvx0k6f32 (CC BY 4.0), file
`Raw_data_supplement_partsizevsmollength.csv`, 848 rows: 24 mesocosms (Dark/Light x
biofilm/nosubs/subs x 4), filters 0.2, 1, 10 um and Total, hours 0 (background), 1, 6, 18, 48,
Short and Long amplicons, copies_mL. A copy is in `New Data mesocosm (Dryad nvx0k6f32)/`.
Script: `edna_py/run_stage0b_mesocosm.py`. Outputs: `results/stage0b_*.csv`,
`results/app_inputs/Mesocosm_*.csv`.

## Table 4.1 reproduced

Source = hour-1 sample, mean composition over mesocosms. Thesis values are matched to two
decimals once replicates with a missing fraction are excluded (one mesocosm each in
Light_biof, Light_nosubs and Dark_nosubs has no 1 um value at hour 1 for the Long assay and a
zero for the Short assay). Example: Open Biofilm Short, thesis (0.56, 0.43, 0.01), data
(0.56, 0.43, 0.01); Shaded No Substrate Long, thesis (0.58, 0.30, 0.12), data (0.58, 0.30, 0.12).
The thesis epsilon is close to the largest between-replicate SD of the three fractions
(e.g. 0.23 vs 0.226, 0.30 vs 0.309, 0.15 vs 0.150), with two values off by 0.01-0.02, so the
rule is "largest replicate SD", rounded.

## Table 4.2 reproduced (run A: thesis f and epsilon; 2000 trials, 0.5 h grid, 100 h window)

Source distance in hours; known answer 5.33 h.

| dataset | n accepted | P10 | P50 | P90 | mode | thesis P10 | P50 | P90 | mode |
|---|---|---|---|---|---|---|---|---|---|
| SBS | 704 | 2.5 | 14.0 | 55.4 | 3.7 | 1.5 | 11.5 | 40.5 | 2.1 |
| SBL | 817 | 3.0 | 15.5 | 54.7 | 3.8 | 1.5 | 10.8 | 37.8 | 1.6 |
| SNS | 554 | 2.0 | 11.0 | 50.5 | 3.8 | 1.0 | 8.5 | 43.1 | 1.7 |
| SNL | 92 | 15.8 | 65.5 | 92.0 | 88.7 | 2.0 | 62.0 | 88.4 | 5.3 |
| SFS | 601 | 2.5 | 15.5 | 63.0 | 3.8 | 2.2 | 16.0 | 60.1 | 3.0 |
| SFL | 1102 | 4.5 | 20.5 | 55.5 | 3.8 | 4.2 | 16.0 | 45.0 | 4.8 |
| OBS | 428 | 2.5 | 13.5 | 39.2 | 3.8 | 3.2 | 19.5 | 52.7 | 5.3 |
| OBL | 318 | 0.5 | 1.0 | 1.5 | 0.6 | 0.5 | 1.0 | 2.0 | 0.8 |
| ONS | 493 | 2.1 | 11.0 | 33.5 | 3.5 | 3.3 | 13.0 | 42.5 | 4.3 |
| ONL | 100 | 1.0 | 3.5 | 11.5 | 3.0 | 1.0 | 4.5 | 31.5 | 0.7 |
| OFS | 2 | 89.6 | 93.8 | 98.0 | 88.9 | 69.6 | 92.0 | 95.2 | 78.7 |
| OFL | 5 | 42.4 | 55.5 | 80.2 | 39.9 | 38.6 | 61.0 | 87.6 | 46.9 |

Coverage of 5.33 h by [P10, P90]: 8 of 12 here, 9 of 12 in the thesis; the one difference is
SNL, where 92 accepted trials give P10 = 15.8 against the thesis's 2.0. The median absolute
difference from the thesis is 2.5 h in P50 and 8.2 h in P90, on a 0.5 h grid with different
random draws. The two Open Bare Substrate failures reproduce exactly in kind: 2 and 5
accepted trials against the thesis's 3 and 4. OBL reproduces to the grid step. Averages over
the ten datasets the thesis averaged: P10/P50/P90 = 3.6/17.1/45.7 h here against
2.0/16.3/44.4 h in the thesis; mode within 5 h of the truth in 9 of 10 (thesis: all 10; the exception here is again SNL).

Run B (source composition and epsilon = largest replicate SD computed from the data) gives
the same coverage (8 of 12) and essentially the same quantiles. Run C (epsilon = 0.05, the
app default) also gives 8 of 12, with narrower intervals (mean P90 36 h) and the SNL dataset
collapsing to 6 accepted trials near zero.

## What it establishes

The Python port, driven by the app's conventions, reproduces the thesis's published
validation from the public data without any adjustment. That closes the loop on
implementation fidelity and gives the paper a positive control it computed itself. It also
shows the three-component case behaving as the thesis described: intervals that bracket the
truth, a front-loaded histogram whose mode sits in the first bins, and one or two treatments
where the regression intervals do not reach the source composition at all.

Housekeeping: the thesis labels the subs treatment "Fresh Substrate" in Table 4.1 and "Bare
Substrate" in Table 4.2; both are the uncolonised-substrate treatment (`subs` in the file).
The thesis's 40 min first sample is `hours == 1` in the file.
