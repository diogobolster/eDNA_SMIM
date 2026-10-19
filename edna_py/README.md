# edna_py: Python implementation of the eDNA analysis tool

`original.py` is a function-for-function port of the model code in Daniel Hallack's R Shiny
application (github.com/dhallackla/edna_app, `appv6+.R` and `aux_fun.R`), including the
implementation conventions the thesis text leaves open. `app.py` is a Shiny for Python interface
with the same tabs as the R application, driven by that module.

Requirements: Python 3.10+, numpy, pandas, scipy, mpmath, matplotlib, shiny (and openpyxl to read
the Snyder xlsx files; playwright and pillow only for `make_figure2.py`).

    pip install numpy pandas scipy mpmath matplotlib shiny openpyxl
    python -m shiny run --port 7799 edna_py/app.py   # then open http://127.0.0.1:7799

Scripts (run from the repository root; `run_all.sh` runs them in order):

    python edna_py/run_stage0_template.py app_R results
    python edna_py/run_stage0b_mesocosm.py data/mesocosm_dryad_nvx0k6f32 results
    python edna_py/run_snyder.py data/snyder results
    python edna_py/run_stage3_smim.py data/snyder results
    python edna_py/make_figures.py data/snyder results figures
    python edna_py/make_figure2.py http://127.0.0.1:7799 results figures/fig2_screens   # app must be running

`tests/test_original.py` checks the port against values computed in R (`tests/r_reference.tsv`):
the linear regression with its 95% intervals, the upper incomplete gamma function of negative
order, Eq. 4.15's k, and the r_h inversion, all to better than 1e-6 relative; plus the
conventions of the Monte Carlo (slope clipping, discarded grid-end solutions, seed
reproducibility) and the recomputed thesis validation.
