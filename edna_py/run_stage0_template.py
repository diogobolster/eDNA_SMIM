"""
Stage 0: the port on Daniel's own template (ES&T 2023 suspended-materials mesocosms).
template_sampling_perSize.csv: 'Distance' = 0, 5, 17, 29 hours after the 1 h sample.
template_source.csv: bucket composition. Known answer: 1 h upstream of the first sample.
App defaults: nsim 1000, tmax 100, dt 2, tol 0.05.

Usage: python3 run_stage0_template.py <app_R dir> <out dir>
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

app = sys.argv[1]; out = sys.argv[2]; os.makedirs(out, exist_ok=True)
samp = pd.read_csv(f"{app}/template_sampling_perSize.csv", encoding="utf-8-sig")
src = pd.read_csv(f"{app}/template_source.csv", encoding="utf-8-sig")
d = O.normalise(samp)

NSIM, TMAX, DT, TOL = 2000, 100, 0.5, 0.05   # thesis validation used 0.5 h steps and 1000-2000 trials
rows = []; regs = []
for tr in sorted(d["Treatment"].unique()):
    target = tr.split("_")[-1]
    tab = O.table_var(d, tr, by_size=True)
    tab.insert(0, "Treatment", tr); regs.append(tab)
    f0 = O.source_composition(src[src["Target"] == target])
    mc = O.mc_new(tab, NSIM, DT, TMAX, seed=1)
    sol = O.df_close(mc, f0["mean"].to_numpy(), TOL, TMAX, denorm=O.denorm_factor(d, tr))
    s = O.summarise_solutions(sol)
    # also with tol = SD of the source composition (the thesis's Table 4.1 rule)
    tol_sd = float(np.nanmean(f0["std"].to_numpy()[:-1])) if len(f0) > 1 else np.nan
    sol2 = O.df_close(mc, f0["mean"].to_numpy(), tol_sd, TMAX, denorm=O.denorm_factor(d, tr))
    s2 = O.summarise_solutions(sol2)
    rows.append(dict(Treatment=tr, known_h=1.0, tol=TOL, **{k: v for k, v in s.items() if k.startswith(("n_sol", "x_"))},
                     tol_sd=round(tol_sd, 3), n_sol_tolsd=s2["n_sol"], x_P10_tolsd=s2["x_P10"], x_P50_tolsd=s2["x_P50"],
                     x_P90_tolsd=s2["x_P90"], x_mode_tolsd=s2["x_mode"]))
R = pd.DataFrame(rows)
pd.concat(regs).to_csv(f"{out}/stage0_template_regressions.csv", index=False, float_format="%.5g")
R.to_csv(f"{out}/stage0_template_prediction.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(R.round(3).to_string())
print("\nknown answer 1 h within [P10,P90] (tol 0.05):", int(((R.x_P10 <= 1) & (R.x_P90 >= 1)).sum()), "of", len(R))
print("known answer 1 h within [P10,P90] (tol = source SD):", int(((R.x_P10_tolsd <= 1) & (R.x_P90_tolsd >= 1)).sum()), "of", len(R))
