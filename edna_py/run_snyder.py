"""
Stages 1, 2 and 4a on the Snyder data with the original app's model (original.py).

Stage 1  app-format CSVs in <out>/app_inputs/ (Distance = 0 at the first station,
         as the app and thesis require; non-detects written as 0, the app's convention).
Stage 2  the app's regressions (replicate-max-normalised, pooled replicates) per
         dataset for the two size fractions and for total; LOD distances at the app
         default (0.05) and at the ddPCR floor; Vdep = lambda Q / w for comparison
         with Snyder.
Stage 4a reference-station prediction: source = composition at the first station
         (mean, SD over its replicates), sampling data = the other four stations
         re-zeroed at the second, known answer = the station spacing (10 m LEEF,
         20 m Fort Custer). App defaults nsim 1000, tmax 100, dt 2, tol 0.05, plus
         tol = SD of the first-station fraction (the thesis's Table 4.1 rule).

Usage: python3 run_snyder.py <data/snyder dir> <out dir>
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

newdata, out = sys.argv[1], sys.argv[2]
os.makedirs(f"{out}/app_inputs", exist_ok=True)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)

LEEF = f"{newdata}/LEEF species specific eDNA removal eDNA concentrations and metadata.xlsx"
FC = f"{newdata}/Fort Custer seasonal eDNA concentration data & metadata.xlsx"

# ---------------------------------------------------------------- Stage 1
def leef_frames():
    tot = pd.read_excel(LEEF, sheet_name="LEEF2021_BothSpeciesTotal", na_values=["NA"])
    c = pd.read_excel(LEEF, sheet_name="LEEF2021_PltData_CARPSize_ForBa", na_values=["NA"])
    s = pd.read_excel(LEEF, sheet_name="LEEF2021_PltData_STLHDSize_ForB", na_values=["NA"])
    bs = pd.concat([c, s], ignore_index=True)
    bs["Treatment"] = bs["Phase"] + "_" + bs["Substrate"]
    size = pd.DataFrame(dict(Treatment=bs["Treatment"], Replicate=bs["Rep"], Filter=bs["Size"].astype(str),
                             Target=bs["Target"], Concentration=bs["GC.mL"].fillna(0.0),
                             Distance=bs["Distance"] - 10, Station_m=bs["Distance"], Q=bs["Q"], w=bs["w"]))
    tot["Treatment"] = tot["Phase"] + "_" + tot["Substrate"]
    total = pd.DataFrame(dict(Treatment=tot["Treatment"], Replicate=tot["Rep"], Target=tot["Target"],
                              Concentration=tot["GC.mL"].fillna(0.0), Distance=tot["Distance"] - 10,
                              Station_m=tot["Distance"], Q=tot["Q"], w=tot["w"]))
    return size, total, 10


def fc_frames():
    ft = pd.read_excel(FC, sheet_name="eDNApltdata_total", na_values=["NA"])
    fs = pd.read_excel(FC, sheet_name="eDNApltdata_bysize", na_values=["NA"]).drop_duplicates()
    size = pd.DataFrame(dict(Treatment=fs["date"], Replicate="R" + fs["rep"].astype(str),
                             Filter=fs["size"].map({1: "1.2", 0: "0.4"}), Target="STHD",
                             Concentration=fs["plteDNA"].fillna(0.0), Distance=fs["station"] - 20,
                             Station_m=fs["station"], Q=fs["Q"], w=fs["w"]))
    total = pd.DataFrame(dict(Treatment=ft["date"], Replicate="R" + ft["rep"].astype(str), Target="STHD",
                              Concentration=ft["plteDNA"].fillna(0.0), Distance=ft["station"] - 20,
                              Station_m=ft["station"], Q=ft["Q"], w=ft["w"]))
    return size, total, 20


EXPS = {"LEEF2021": leef_frames(), "FortCuster": fc_frames()}
for name, (size, total, _) in EXPS.items():
    size.drop(columns=["Station_m", "Q", "w"]).to_csv(f"{out}/app_inputs/{name}_sampling_perSize.csv", index=False)
    total.drop(columns=["Station_m", "Q", "w"]).to_csv(f"{out}/app_inputs/{name}_sampling_Total.csv", index=False)

# ---------------------------------------------------------------- Stage 2
FLOOR = {"CARP": 0.0351, "STLHD": 0.0557, "STHD": 0.0557}   # smallest quantified GC/mL in the files
reg_rows, lod_rows = [], []
for name, (size, total, spacing) in EXPS.items():
    for target in size["Target"].unique():
        ds = O.normalise(size[size["Target"] == target])
        dt_ = O.normalise(total[total["Target"] == target])
        for tr in sorted(ds["Treatment"].unique()):
            Q = size.loc[size["Treatment"] == tr, "Q"].iloc[0]; w = size.loc[size["Treatment"] == tr, "w"].iloc[0]
            tabs = [O.table_var(ds, tr, by_size=True), O.table_var(dt_, tr, by_size=False)]
            for tab, dd in zip(tabs, [ds, dt_]):
                if len(tab) == 0:
                    continue
                for f in tab["Filter"].unique():
                    ic = tab[(tab.Filter == f) & (tab.term == "(Intercept)")].iloc[0]
                    sl = tab[(tab.Filter == f) & (tab.term == "Distance")].iloc[0]
                    lam = -sl["estimate"]
                    reg_rows.append(dict(experiment=name, target=target, dataset=tr, component=f, n=int(sl["n"]),
                                         slope=sl["estimate"], slope_low=sl["conf.low"], slope_high=sl["conf.high"],
                                         lambda_=lam, lambda_CI_excludes_0=sl["conf.high"] < 0,
                                         intercept=ic["estimate"], r2=sl["r.squared"], p_slope=sl["p.value"],
                                         Sw_m=(1 / lam if lam > 0 else np.inf), Vdep_m_per_day=lam * Q / w,
                                         denorm_factor=O.denorm_factor(dd, tr)))
                mf = O.max_factor(dd, tr)
                for lod_name, lod in [("app_default_0.05", 0.05), ("ddPCR_floor", FLOOR[target])]:
                    L = O.lod_distance(tab, mf, lod)
                    for _, r in L.iterrows():
                        lod_rows.append(dict(experiment=name, target=target, dataset=tr, component=r["Filter"],
                                             LOD=lod_name, LOD_value=lod, x_LOD_avg=r["avg"], x_LOD_upper=r["upper"],
                                             x_LOD_lower=r["lower"], first_station_m=spacing))
REG = pd.DataFrame(reg_rows).rename(columns={"lambda_": "lambda"})
LODT = pd.DataFrame(lod_rows)
for c in ["x_LOD_avg", "x_LOD_upper", "x_LOD_lower"]:
    LODT[c + "_from_injection"] = LODT[c] + LODT["first_station_m"]
REG.to_csv(f"{out}/stage2_regressions_original_method.csv", index=False, float_format="%.5g")
LODT.to_csv(f"{out}/stage2_lod_distances.csv", index=False, float_format="%.5g")

# ---------------------------------------------------------------- Stage 4a
NSIM, TMAX, DT, TOL = 1000, 100, 2, 0.05
pred_rows = []
for name, (size, total, spacing) in EXPS.items():
    for target in size["Target"].unique():
        sub = size[size["Target"] == target]
        for tr in sorted(sub["Treatment"].unique()):
            g = sub[sub["Treatment"] == tr]
            # source = first station; sampling = the rest, re-zeroed at the second station
            src = g[g["Station_m"] == g["Station_m"].min()]
            src = src[src["Concentration"] > 0]
            # as f_1st_sample() does (frac < 1): a replicate with only one fraction quantified is dropped
            nf_rep = src.groupby("Replicate")["Filter"].transform("nunique")
            src = src[nf_rep == g["Filter"].nunique()]
            samp = g[g["Station_m"] > g["Station_m"].min()].copy()
            samp["Distance"] = samp["Station_m"] - (g["Station_m"].min() + spacing)
            d = O.normalise(samp)
            tab = O.table_var(d, tr, by_size=True)
            f0 = O.source_composition(src)                     # mean/SD over the 5 reps at station 1
            if tab["Filter"].nunique() < 2 or len(f0) < 2:
                pred_rows.append(dict(experiment=name, target=target, dataset=tr, note="fewer than 2 fitted components"))
                continue
            mc = O.mc_new(tab, NSIM, DT, TMAX, seed=1)
            dn = O.denorm_factor(d, tr)
            tol_sd = float(f0["std"].iloc[0])                  # SD of the first (n-1 = 1) fraction
            for tol_name, tol in [("app_0.05", TOL), ("source_SD", tol_sd)]:
                sol = O.df_close(mc, f0["mean"].to_numpy(), tol, TMAX, denorm=dn)
                s = O.summarise_solutions(sol)
                known = spacing
                pred_rows.append(dict(experiment=name, target=target, dataset=tr, tol_rule=tol_name, tol=tol,
                                      f_source_first=f0["mean"].iloc[0], f_source_sd=tol_sd,
                                      known_x_m=known, **s,
                                      covered=(s["x_P10"] <= known <= s["x_P90"]) if s["n_sol"] > 1 else False,
                                      rel_err_mode=(s["x_mode"] - known) / known if s["n_sol"] > 1 else np.nan,
                                      rel_width=(s["x_P90"] - s["x_P10"]) / TMAX if s["n_sol"] > 1 else np.nan,
                                      pred_source_from_injection_mode=(2 * spacing - s["x_mode"]) if s["n_sol"] > 1 else np.nan))
PRED = pd.DataFrame(pred_rows)
PRED.to_csv(f"{out}/stage4a_reference_station_prediction.csv", index=False, float_format="%.4g")

# ---------------------------------------------------------------- summary to stdout
print("STAGE 2: lambda 95% CI excludes 0 (original method, normalised)")
print(REG.groupby(["experiment", "component"])["lambda_CI_excludes_0"].agg(["sum", "count"]))
print("\nFort Custer, original method:")
print(REG[REG.experiment == "FortCuster"][["dataset", "component", "n", "lambda", "r2", "Sw_m", "Vdep_m_per_day", "lambda_CI_excludes_0"]].round(4).to_string())
print("\nSTAGE 4a summary by tolerance rule:")
ok = PRED[PRED.get("tol_rule").notna()] if "tol_rule" in PRED else PRED
for (e, t), g in ok.groupby(["experiment", "tol_rule"]):
    g2 = g[g.n_sol > 1]
    print(f"{e:10s} {t:10s} datasets {len(g)}, with >1 solution {len(g2)}, truth coverage {g2.covered.mean() if len(g2) else float('nan'):.2f}, "
          f"median |rel err mode| {g2.rel_err_mode.abs().median() if len(g2) else float('nan'):.2f}, median rel width {g2.rel_width.median() if len(g2) else float('nan'):.2f}, "
          f"median n_sol {g2.n_sol.median() if len(g2) else float('nan'):.0f}")
print("\nFort Custer 4a detail:")
print(ok[ok.experiment == "FortCuster"][["dataset", "tol_rule", "tol", "f_source_first", "n_sol", "x_P10", "x_P50", "x_P90", "x_mode", "known_x_m", "covered"]].round(3).to_string())
