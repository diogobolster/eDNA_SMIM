"""
Step 0 baseline: a plain semi-log regression (no replicate-max normalisation; superseded by edna_py/original.py) applied
to every Snyder dataset, nothing more.

For each dataset and component: lm(log C ~ x) on all replicate points with C > 0,
lambda = -slope, C0 = exp(intercept), 95% CI on both, lambda truncated at 0 exactly
as the app does. Non-detects (NA or 0) are dropped, exactly as the app does.

Usage:  python3 baseline_fits.py "<path to New Data folder>"  [output.csv]
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

newdata = sys.argv[1] if len(sys.argv) > 1 else "New Data"
out = sys.argv[2] if len(sys.argv) > 2 else "03_baseline_fits.csv"

LEEF = f"{newdata}/LEEF species specific eDNA removal eDNA concentrations and metadata.xlsx"
FC = f"{newdata}/Fort Custer seasonal eDNA concentration data & metadata.xlsx"


def fit_one(x, conc, conf=0.95):
    """Same as fit_one() in edna-app/R/model.R."""
    x = np.asarray(x, float); conc = np.asarray(conc, float)
    keep = np.isfinite(x) & np.isfinite(conc) & (conc > 0)
    x, conc = x[keep], conc[keep]
    if len(x) < 3 or len(np.unique(x)) < 2:
        return None
    r = stats.linregress(x, np.log(conc))
    n = len(x)
    tcrit = stats.t.ppf(0.5 + conf / 2, n - 2)
    slope_lo, slope_up = r.slope - tcrit * r.stderr, r.slope + tcrit * r.stderr
    icept_lo, icept_up = r.intercept - tcrit * r.intercept_stderr, r.intercept + tcrit * r.intercept_stderr
    return dict(
        n=n,
        lambda_=max(-r.slope, 0),
        lambda_low=max(-slope_up, 0),
        lambda_up=max(-slope_lo, 0),
        slope_raw=r.slope,                 # untruncated, so flat/rising profiles are visible
        C0=np.exp(r.intercept), C0_low=np.exp(icept_lo), C0_up=np.exp(icept_up),
        r2=r.rvalue ** 2, p_slope=r.pvalue,
    )


rows = []

# ---- LEEF 2021 (Snyder Ch. 5): 4 streams x 4 phases x 2 species, by size and total
tot = pd.read_excel(LEEF, sheet_name="LEEF2021_BothSpeciesTotal", na_values=["NA"])
carp = pd.read_excel(LEEF, sheet_name="LEEF2021_PltData_CARPSize_ForBa", na_values=["NA"])
sthd = pd.read_excel(LEEF, sheet_name="LEEF2021_PltData_STLHDSize_ForB", na_values=["NA"])
bysize = pd.concat([carp, sthd], ignore_index=True)

for (tg, ph, sub), g in bysize.groupby(["Target", "Phase", "Substrate"]):
    for sz, gg in g.groupby("Size"):
        f = fit_one(gg["Distance"], gg["GC.mL"])
        if f:
            rows.append(dict(experiment="LEEF2021", target=tg, dataset=f"{ph}_{sub}", component=f"{sz}um",
                             n_nondetect=int(gg["GC.mL"].isna().sum()), Q=gg["Q"].iloc[0], w=gg["w"].iloc[0], **f))
for (tg, ph, sub), g in tot.groupby(["Target", "Phase", "Substrate"]):
    f = fit_one(g["Distance"], g["GC.mL"])
    if f:
        rows.append(dict(experiment="LEEF2021", target=tg, dataset=f"{ph}_{sub}", component="total",
                         n_nondetect=int((g["GC.mL"] == 0).sum()), Q=g["Q"].iloc[0], w=g["w"].iloc[0], **f))

# ---- Fort Custer (Snyder Ch. 6): 4 seasons, Steelhead only, conductivity-corrected
ft = pd.read_excel(FC, sheet_name="eDNApltdata_total", na_values=["NA"])
fs = pd.read_excel(FC, sheet_name="eDNApltdata_bysize", na_values=["NA"]).drop_duplicates()
for dt, g in ft.groupby("date"):
    f = fit_one(g["station"], g["plteDNA"])
    rows.append(dict(experiment="FortCuster", target="STHD", dataset=dt, component="total",
                     n_nondetect=int(g["plteDNA"].isna().sum()), Q=g["Q"].iloc[0], w=g["w"].iloc[0], **f))
for (dt, sz), g in fs.groupby(["date", "size"]):
    f = fit_one(g["station"], g["plteDNA"])
    rows.append(dict(experiment="FortCuster", target="STHD", dataset=dt, component={1: "1.2um", 0: "0.4um"}[sz],
                     n_nondetect=int(g["plteDNA"].isna().sum()), Q=g["Q"].iloc[0], w=g["w"].iloc[0], **f))

R = pd.DataFrame(rows).rename(columns={"lambda_": "lambda"})
R["Sw_m"] = np.where(R["lambda"] > 0, 1 / R["lambda"], np.inf)
R["lambda_CI_excludes_0"] = R["lambda_low"] > 0
R["Vdep_m_per_day"] = R["lambda"] * R["Q"] / R["w"]     # Q is m3/day in both files
R.to_csv(out, index=False, float_format="%.5g")

print(R.to_string())
print("\nFits with lambda 95% CI excluding 0:")
print(R.groupby(["experiment", "component"])["lambda_CI_excludes_0"].agg(["sum", "count"]))
