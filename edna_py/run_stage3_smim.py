"""
Stage 3: SMIM-rates tab (Eq. 4.15, Sec. 4.4) on LEEF 2021, the original app's way.

Inputs
  k_data  : -slope of the app's regression (Stage 2 output), per stream x phase x species
            x {1.2, 0.4, total}, in 1/m, with its 95% CI.
  U, D, Lambda, beta : rhodamine breakthrough fits per stream x phase from the covariables
            file (Snyder Ch. 5 Table / 'LEEF2021Ancillary').
  t1, t2  : hard-coded in appv6+.R, 1 s and 1e5 s.
  r_s     : water-column removal rate per species and size, fitted here from the Dryad
            polydisperse mesocosm time series (Brandao-Dias et al. 2023a, treatment AF =
            sequential filtration of unfiltered source water, 1/6/18/30 h), semi-log
            regression pooling replicates, converted from 1/h to 1/s. Size mapping:
            LEEF '1.2' (>1.2 um) ~ Dryad 10 um + 1 um filters summed; LEEF '0.4'
            (0.4-1.2 um) ~ Dryad 0.4 um filter; 'total' ~ treatment N (no filtration).
Output
  r_h (1/s) per combination by least squares on |k_model - k_data|, plus the r_h implied
  by the two ends of the k CI, the k the model gives with r_h = 0 (water column only), and
  sensitivity to t2 (1e3, 1e5) and to r_s (0, fitted, 2x fitted).

Usage: python3 run_stage3_smim.py <data/snyder dir> <out dir>   (needs <out>/stage2_regressions_original_method.csv)
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

newdata, out = sys.argv[1], sys.argv[2]
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)

# ---------------------------------------------------------------- r_s from the Dryad mesocosm series
dry = pd.read_csv(f"{newdata}/doi_10_5061_dryad_37pvmcvnw__v20221219/polidisperse.experiment.datasheet.csv")
dry = dry[dry.hours > 0]


def fit_rate(x, y):
    """semi-log regression of C on t (pooled replicates, zeros dropped): rate = -slope (1/h) with 95% CI"""
    y = np.asarray(y, float); x = np.asarray(x, float)
    k = y > 0; x, y = x[k], np.log(y[k])
    r = stats.linregress(x, y); t = stats.t.ppf(0.975, len(x) - 2)
    return -r.slope, -(r.slope + t * r.stderr), -(r.slope - t * r.stderr), r.rvalue ** 2, len(x)


rs_rows = []
for tg in ["CARP", "STHD"]:
    af = dry[(dry.target == tg) & (dry.treatment == "AF")]
    big = af[af["filter"].isin([10.0, 1.0])].groupby(["replicate", "hours"])["copies.mL"].sum().reset_index()
    small = af[af["filter"] == 0.4]
    tot = dry[(dry.target == tg) & (dry.treatment == "N")]
    for comp, g in [("1.2", big), ("0.4", small), ("total", tot)]:
        rate, lo, hi, r2, n = fit_rate(g.hours, g["copies.mL"])
        rs_rows.append(dict(target={"STHD": "STLHD"}.get(tg, tg), component=comp, rs_per_h=rate, rs_low_per_h=lo, rs_high_per_h=hi,
                            r2=r2, n=n, rs_per_s=rate / 3600, source="Dryad AF (10+1 um)" if comp == "1.2" else ("Dryad AF 0.4 um" if comp == "0.4" else "Dryad N (unfiltered)")))
RS = pd.DataFrame(rs_rows)
RS.to_csv(f"{out}/stage3_rs_from_mesocosm.csv", index=False, float_format="%.4g")
print("r_s from the Dryad mesocosm series (1/h):\n", RS.round(4).to_string(), "\n")

# ---------------------------------------------------------------- conservative parameters
cov = pd.read_excel(f"{newdata}/LEEF species specific eDNA removal environmental covariables and metadata.xlsx", sheet_name="LEEF2021Ancillary")
phase_map = {("early", "shade"): "S", ("early", "unshade"): "LB", ("late", "unshade"): "HB", ("late", "shade"): "RS"}
sub_map = {"Cobble": "COB", "Pea Gravel": "PG", "Sand": "SAND", "Mix": "MIX"}
cov["dataset"] = [phase_map[(b, s)] + "_" + sub_map[sub] for b, s, sub in zip(cov["biofilm.text"], cov["shade.text"], cov["substrate"])]
cov = cov.set_index("dataset")[["U", "D", "Lambda", "Beta", "U_SE", "D_SE", "Lambda_SE", "Beta_SE"]]

# ---------------------------------------------------------------- inversion
REG = pd.read_csv(f"{out}/stage2_regressions_original_method.csv")
REG = REG[REG.experiment == "LEEF2021"].copy()
REG["component"] = REG["component"].astype(str)
rows = []
for _, r in REG.iterrows():
    c = cov.loc[r.dataset]
    rs_row = RS[(RS.target == r.target) & (RS.component == r.component)].iloc[0]
    rs = rs_row.rs_per_s
    k_data, k_lo, k_hi = -r.slope, -r.slope_high, -r.slope_low        # 1/m; k_lo <= k_data <= k_hi
    U, D, L, B = c.U, c.D, c.Lambda, c.Beta
    k_rh0 = O.k_func(rs, 0.0, U, D, L, B)
    def inv(k, rs_=rs, logT2=O.LOGT2):
        if not np.isfinite(k) or k <= 0:
            return np.nan, np.nan
        return O.fit_rh(k, rs_, U, D, L, B, logT2=logT2)
    rh, res = inv(k_data)
    rh_lo, _ = inv(k_lo); rh_hi, _ = inv(k_hi)
    rh_t2e3, _ = inv(k_data, logT2=3.0)
    rh_rs0, _ = inv(k_data, rs_=0.0)
    rh_rs2, _ = inv(k_data, rs_=2 * rs)
    rows.append(dict(target=r.target, dataset=r.dataset, component=r.component, n=r.n, r2=r.r2,
                     k_data_per_m=k_data, k_CI_low=k_lo, k_CI_high=k_hi, k_CI_excludes_0=r.lambda_CI_excludes_0,
                     U=U, D=D, Lambda=L, beta=B, rs_per_s=rs,
                     k_if_rh_0=k_rh0, water_column_explains_k=(k_rh0 >= k_data),
                     rh_per_s=rh, rh_residual=res, rh_at_k_CI_low=rh_lo, rh_at_k_CI_high=rh_hi,
                     rh_over_rs=(rh / rs if np.isfinite(rh) and rs > 0 else np.nan),
                     rh_t2_1e3=rh_t2e3, rh_rs_0=rh_rs0, rh_rs_2x=rh_rs2))
SM = pd.DataFrame(rows)
SM.to_csv(f"{out}/stage3_smim_rh.csv", index=False, float_format="%.4g")

ok = SM[SM.k_CI_excludes_0 & np.isfinite(SM.rh_per_s)]
print("Inversions attempted:", len(SM), " k > 0:", int((SM.k_data_per_m > 0).sum()), " slope CI excludes 0:", int(SM.k_CI_excludes_0.sum()))
print("Water column (r_s alone) explains k in", int(SM.water_column_explains_k.sum()), "of", len(SM), "cases")
print("\nr_h (1/s) where the slope CI excludes 0:")
print(ok[["target", "dataset", "component", "k_data_per_m", "rs_per_s", "k_if_rh_0", "rh_per_s", "rh_at_k_CI_low", "rh_at_k_CI_high", "rh_over_rs", "rh_t2_1e3", "rh_rs_0"]].to_string(float_format=lambda v: f"{v:.3g}"))
print("\nmedian r_h by component (significant fits):")
print(ok.groupby(["target", "component"])["rh_per_s"].agg(["median", "min", "max", "count"]).to_string(float_format=lambda v: f"{v:.3g}"))
print("\nmedian r_h / r_s:", ok.rh_over_rs.median())
print("sensitivity: median rh(t2=1e3)/rh(t2=1e5) =", (ok.rh_t2_1e3 / ok.rh_per_s).median(), "; median rh(rs=0)/rh =", (ok.rh_rs_0 / ok.rh_per_s).median())
