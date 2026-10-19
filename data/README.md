# Data

All files are included so that every number in the paper can be regenerated with `run_all.sh`.
Each dataset keeps the licence of its source.

## `snyder/`: field additions in flowing water (Snyder et al.)

Steady-state eDNA additions with a known source, sampled downstream by size fraction
(1.2 um and 0.4 um filters, sequential) and in total.

| File | Content |
|---|---|
| `LEEF species specific eDNA removal eDNA concentrations and metadata.xlsx` | ND-LEEF 2021: four substrate streams, four phases (early/late open, early/late closed), common carp and steelhead eDNA, stations at 10 to 50 m |
| `LEEF species specific eDNA removal environmental covariables and metadata.xlsx` | Velocity U, dispersion D, SMIM Lambda and beta per stream and phase, from conservative tracer tests |
| `Fort Custer seasonal eDNA concentration data & metadata.xlsx` | Fort Custer stream: four seasons, stations at 20 to 100 m, conductivity-corrected concentrations |
| `Fort Custer seasonal environmental covariables data & metadata.xlsx` | Covariates per season |

Source: Snyder, E. D. (2024). *Using Experimental and Field Approaches to Identify Predictors and
Drivers of eDNA Transport and Removal in Streams and Rivers*. PhD thesis, University of Notre Dame,
https://doi.org/10.7274/27071209.v1, and the associated papers cited in the manuscript. Included
with the permission of the authors; please cite the thesis or papers when reusing these data.

`doi_10_5061_dryad_37pvmcvnw__v20221219/polidisperse.experiment.datasheet.csv`: the polydisperse
mesocosm time series used to set the surface-water decay rate r_s for the SMIM inversion
(Brandão-Dias et al., Dryad https://doi.org/10.5061/dryad.37pvmcvnw).

## `mesocosm_dryad_nvx0k6f32/`: recirculating mesocosm experiment

`Raw_data_supplement_partsizevsmollength.csv`, `Metadata.csv` and the Dryad `README.md`, from
Brandão-Dias et al. (2025), *Comparing the fate of eDNA by particle sizes and molecule lengths in
recirculating streams*, Dryad https://doi.org/10.5061/dryad.nvx0k6f32 (CC BY 4.0). This is the
experiment behind the source-localisation validation in Hallack (2025), Tables 4.1 and 4.2, which
`edna_py/run_stage0b_mesocosm.py` recomputes.

## `../app_R/template_*.csv`

The example inputs shipped with the R application (ES&T 2023 suspended-materials mesocosms:
Brandão-Dias et al. 2023, *Environ. Sci. Technol.* 57, 13161).

## `../results/app_inputs/`

The same datasets rewritten in the application's input format, one file per dataset and
fraction type, so that any analysis in the paper can be repeated by loading the file in either
application.
