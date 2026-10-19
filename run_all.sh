#!/usr/bin/env bash
# Reproduces results/ and figures/ from data/ with the published method, in order.
set -euo pipefail
cd "$(dirname "$0")"
python edna_py/run_stage0_template.py app_R results
python edna_py/run_stage0b_mesocosm.py data/mesocosm_dryad_nvx0k6f32 results
python edna_py/run_snyder.py data/snyder results
python edna_py/run_stage3_smim.py data/snyder results
python edna_py/make_figures.py data/snyder results figures
echo "Figure 2 needs the app running: python -m shiny run --port 7799 edna_py/app.py &"
echo "then: python edna_py/make_figure2.py http://127.0.0.1:7799 results figures/fig2_screens"
