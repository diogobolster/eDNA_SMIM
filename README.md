# eDNA_SMIM: interpreting environmental DNA concentration profiles in streams

A web application and model code for four analyses of eDNA concentration profiles in flowing
water: removal rates per component (size class), limit-of-detection distance, hyporheic and
surface-water reaction rates (reactive stochastic mobile-immobile model, SMIM), and source
distance from the change in particle-size composition downstream. The models are those of
Hallack (2025), Chapter 4. Two implementations are provided:

| Folder | What | Run |
|---|---|---|
| `app_R/` | The original R Shiny application (Daniel Hallack, 2024; copied from github.com/dhallackla/edna_app at commit 4ac13a7) | `library(shiny); runApp("app_R")` |
| `edna_py/` | A Python port of the model routines, function for function (`original.py`), and a Shiny for Python application with the same tabs (`app.py`) | `python -m shiny run --port 7799 edna_py/app.py` |

The two implementations agree on every model routine (see `tests/`). Everything else in the
repository reproduces the evaluation reported in the accompanying paper (Hallack et al., in
preparation for *Environmental Modelling & Software*): the data, the scripts that run the
published method on them, the result tables and the figures.

## Install (Python)

Python 3.10 or later.

    pip install -r requirements.txt
    python -m pytest tests          # checks the port against the R reference values

## Run the application

    python -m shiny run --port 7799 edna_py/app.py
    # open http://127.0.0.1:7799

Input CSVs follow the R application's templates in `app_R/` (`template_sampling_perSize.csv`,
`template_sampling_Total.csv`, `template_source.csv`): columns `Treatment, Replicate, Filter,
Target, Concentration, Distance` (or `Hours`). Ready-made inputs for the field and mesocosm
datasets are in `results/app_inputs/`. The R application needs the packages listed in
`app_R/README.md`.

## Reproduce the paper

`run_all.sh` runs every stage in order and rewrites `results/` and `figures/`:

    bash run_all.sh

The stages, each a script in `edna_py/`:

| Stage | Script | Input | Output |
|---|---|---|---|
| 0 | `run_stage0_template.py` | `app_R/template_*.csv` | `results/stage0_template_*.csv` |
| 0b | `run_stage0b_mesocosm.py` | `data/mesocosm_dryad_nvx0k6f32/` | `results/stage0b_*.csv`, `results/app_inputs/Mesocosm_*.csv` (thesis Tables 4.1 and 4.2 recomputed) |
| 1, 2, 4a | `run_snyder.py` | `data/snyder/` | regressions, LOD distances, reference-station source predictions, `results/app_inputs/{LEEF2021,FortCuster}_*.csv` |
| 3 | `run_stage3_smim.py` | `data/snyder/`, stage 2 output | `results/stage3_*.csv` (SMIM inversion for r_h) |
| figures | `make_figures.py` | `results/` | `figures/fig1, fig3 to fig6` |
| figure 2 | `make_figure2.py` | a running app | `figures/fig2_screens/`, `figures/fig2_app_screens.png` (needs playwright and a Chromium) |

The Monte Carlo stages use fixed seeds; quantiles move by about one grid step between seeds and
between the R and Python random number generators.

## Data

See `data/README.md` for provenance and licences. In short: the ND-LEEF 2021 and Fort Custer
field datasets (Snyder et al.), the Brandão-Dias et al. (2025) mesocosm dataset (Dryad
doi:10.5061/dryad.nvx0k6f32, CC BY 4.0) and the polydisperse mesocosm series (Dryad
doi:10.5061/dryad.37pvmcvnw).

## Documentation

`docs/` holds the working notes written while the method was applied as published: what the
thesis and the code state (`01`), the data inventory (`02`), the plan (`04`), where the code
departs from the thesis text (`06`), and the results of each stage (`07`, `08`, `10`).

## Citation

See `CITATION.cff`. The models are described in: Hallack, D. M. C. (2025). *Mixing and Removal
Processes in Complex Environmental Flow Systems* (PhD thesis, University of Notre Dame),
Chapter 4. https://doi.org/10.7274/28756664.v1

## Licence

MIT, see `LICENSE`. The data files keep their own licences (see `data/README.md`).
