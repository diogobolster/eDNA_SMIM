"""
Stage 0b: recompute the thesis validation (Hallack 2025, Sec. 4.3.1.2, Tables 4.1-4.2) with
the Python port, from the archived mesocosm data (Brandao-Dias et al. 2025; Dryad
doi:10.5061/dryad.nvx0k6f32, file Raw_data_supplement_partsizevsmollength.csv).

Design exactly as the thesis states it:
  source   = the first post-release sample (hours == 1 in the file; 40 min in the thesis),
             composition f^s = mean over the four mesocosms of each treatment, replicates with
             a missing fraction excluded (this reproduces Table 4.1 to two decimals);
  epsilon  = thesis Table 4.1 value (run A) and, as a check, the largest between-replicate
             SD of the three fractions (run B);
  sampling = hours 6, 18, 48 shifted to 0, 12, 42; known answer 5.33 h;
  components = the three filter sizes 0.2, 1, 10 um (thesis case C);
  settings = x_max 100 h, step 0.5 h, 2000 trials; app conventions throughout.

Usage: python3 run_stage0b_mesocosm.py "<folder with the Dryad csv>" <out dir>
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

src_dir, out = sys.argv[1], sys.argv[2]
os.makedirs(f"{out}/app_inputs", exist_ok=True)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)

d = pd.read_csv(f"{src_dir}/Raw_data_supplement_partsizevsmollength.csv")
d = d[(~d["Background"]) & (d["filter"] != "Total") & (d["treatment"] != "Bucket")].copy()
d["Treatment"] = d["treatment"] + "_" + d["Target"]

# thesis Table 4.1: acronym, f1 (0.2), f2 (1), f3 (10), epsilon; Ans = 5.33 h
T41 = {
    "Dark_biof_Short": ("SBS", 0.65, 0.34, 0.01, 0.08), "Dark_biof_Long": ("SBL", 0.73, 0.26, 0.01, 0.11),
    "Dark_nosubs_Short": ("SNS", 0.58, 0.32, 0.10, 0.23), "Dark_nosubs_Long": ("SNL", 0.58, 0.30, 0.12, 0.30),
    "Dark_subs_Short": ("SFS", 0.53, 0.44, 0.03, 0.15), "Dark_subs_Long": ("SFL", 0.58, 0.40, 0.02, 0.25),
    "Light_biof_Short": ("OBS", 0.56, 0.43, 0.01, 0.02), "Light_biof_Long": ("OBL", 0.78, 0.21, 0.01, 0.02),
    "Light_nosubs_Short": ("ONS", 0.55, 0.34, 0.11, 0.09), "Light_nosubs_Long": ("ONL", 0.53, 0.38, 0.09, 0.14),
    "Light_subs_Short": ("OFS", 0.59, 0.39, 0.02, 0.06), "Light_subs_Long": ("OFL", 0.66, 0.31, 0.02, 0.18),
}
# thesis Table 4.2 (source distance, h): P10, P50, P90, mode
T42 = {"SBS": (1.5, 11.5, 40.5, 2.1), "SBL": (1.5, 10.8, 37.8, 1.6), "SNS": (1, 8.5, 43.1, 1.7), "SNL": (2, 62, 88.4, 5.3),
       "SFS": (2.2, 16, 60.1, 3), "SFL": (4.2, 16, 45, 4.8), "OBS": (3.2, 19.5, 52.7, 5.3), "OBL": (0.5, 1, 2, 0.8),
       "ONS": (3.3, 13, 42.5, 4.3), "ONL": (1, 4.5, 31.5, 0.7), "OFS": (69.6, 92, 95.2, 78.7), "OFL": (38.6, 61, 87.6, 46.9)}
# note: thesis Table 4.2 labels the last two rows "Open Bare Substrate"; Fresh/Bare both mean uncolonised substrate (subs).

samp = d[d["hours"].isin([6, 18, 48])].copy()
samp_app = pd.DataFrame(dict(Treatment=samp["Treatment"], Replicate="M" + samp["mesocosm"].astype(str), Filter=samp["filter"],
                             Target=samp["Target"], Concentration=samp["copies_mL"].fillna(0.0), Distance=samp["hours"] - 6))
src = d[d["hours"] == 1].copy()
src_app = pd.DataFrame(dict(Treatment=src["Treatment"], Replicate="M" + src["mesocosm"].astype(str), Filter=src["filter"],
                            Target=src["Target"], Concentration=src["copies_mL"]))
samp_app.to_csv(f"{out}/app_inputs/Mesocosm_sampling_perSize_hours.csv", index=False)
src_app.to_csv(f"{out}/app_inputs/Mesocosm_source.csv", index=False)

KNOWN = 5.33
NSIM, TMAX, DT = 2000, 100, 0.5
dn = O.normalise(samp_app)
rows = []; regs = []
for tr in sorted(dn["Treatment"].unique()):
    acr, f1, f2, f3, eps_thesis = T41[tr]
    tab = O.table_var(dn, tr, by_size=True); tab.insert(0, "Treatment", tr); regs.append(tab)
    # source: replicates with all three fractions quantified
    s = src_app[(src_app["Treatment"] == tr) & (src_app["Concentration"] > 0)]
    nf = s.groupby("Replicate")["Filter"].transform("nunique"); s = s[nf == 3]
    f0 = O.source_composition(s)
    eps_sd = float(f0["std"].max())
    mc = O.mc_new(tab, NSIM, DT, TMAX, seed=1)
    denorm = O.denorm_factor(dn, tr)
    for run, fs, eps in [("A_thesis_f_eps", np.array([f1, f2, f3]), eps_thesis),
                         ("B_data_f_maxSD", f0["mean"].to_numpy(), eps_sd),
                         ("C_data_f_tol0.05", f0["mean"].to_numpy(), 0.05)]:
        sol = O.df_close(mc, fs, eps, TMAX, denorm=denorm)
        st = O.summarise_solutions(sol)
        rows.append(dict(acronym=acr, treatment=tr, run=run, f1=fs[0], f2=fs[1], f3=fs[2], eps=eps, n_src_reps=s["Replicate"].nunique(),
                         known_h=KNOWN, **{k: v for k, v in st.items() if k.startswith(("n_sol", "x_"))},
                         covered=(st["x_P10"] <= KNOWN <= st["x_P90"]) if st["n_sol"] > 1 else False,
                         thesis_P10=T42[acr][0], thesis_P50=T42[acr][1], thesis_P90=T42[acr][2], thesis_mode=T42[acr][3],
                         thesis_covered=(T42[acr][0] <= KNOWN <= T42[acr][2])))
R = pd.DataFrame(rows)
pd.concat(regs).to_csv(f"{out}/stage0b_mesocosm_regressions.csv", index=False, float_format="%.5g")
R.to_csv(f"{out}/stage0b_mesocosm_table42_recomputed.csv", index=False, float_format="%.4g")
for run, g in R.groupby("run"):
    print(f"\n=== {run}: coverage {int(g.covered.sum())}/12 (thesis {int(g.thesis_covered.sum())}/12), "
          f"median n_sol {g.n_sol.median():.0f}, mean P10/P50/P90/mode {g.x_P10.mean():.1f}/{g.x_P50.mean():.1f}/{g.x_P90.mean():.1f}/{g.x_mode.mean():.1f} "
          f"(thesis 2.0/16.3/44.4/2.9)")
    print(g[["acronym", "f1", "f2", "f3", "eps", "n_sol", "x_P10", "x_P50", "x_P90", "x_mode", "covered", "thesis_P10", "thesis_P50", "thesis_P90", "thesis_mode"]].round(2).to_string(index=False))
