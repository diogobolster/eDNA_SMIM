"""
Figures 1, 3, 4, 5, 6 for the EMS manuscript, from the result tables and the model module.
Figure 2 (application screenshots) is made from the running Shiny app, not here.

Usage: python3 make_figures.py <data/snyder dir> <results dir> <figures dir>
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

newdata, res, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)

# validated categorical palette (dataviz reference instance, slots 1-3) + ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 1.6, "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans",
})
COMP_COL = {"1.2": BLUE, "0.4": ORANGE, "total": AQUA}
COMP_LAB = {"1.2": "> 1.2 um", "0.4": "0.4 to 1.2 um", "total": "total"}


def save(fig, name):
    fig.savefig(f"{out}/{name}.png"); fig.savefig(f"{out}/{name}.pdf"); plt.close(fig)
    print("wrote", name)


def panel_label(ax, s):
    ax.text(-0.08, 1.04, s, transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="right")


# ============================================================== Figure 1: concept
def fig1():
    rng = np.random.default_rng(3)
    lam = np.array([0.045, 0.020, 0.004])          # 1/m, components A, B, C
    C0 = np.array([0.30, 0.40, 0.30])              # composition at the source, x = -L
    L = 25.0
    x = np.linspace(-L, 60, 300)
    C = C0[:, None] * np.exp(-lam[:, None] * (x[None, :] + L))
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3), gridspec_kw={"width_ratios": [1.15, 1.2, 0.9]})

    ax = axs[0]
    for i, (c, col, lab) in enumerate(zip(C, [BLUE, ORANGE, AQUA], ["A", "B", "C"])):
        ax.plot(x, c, color=col); ax.text(61, c[-1], lab, color=col, va="center", fontweight="bold")
    stations = np.array([0, 10, 20, 30, 40])
    for s in stations:
        ax.axvline(s, color=INK2, lw=0.6, ls=":")
    ax.axvline(-L, color=INK, lw=1.0); ax.text(-L + 1.5, 0.47, "source\n(composition known)", ha="left", va="top", fontsize=7.5)
    ax.text(20, 0.47, "sampling stations", ha="center", va="top", fontsize=7.5, color=INK2)
    ax.set_xlabel("distance x (m), x = 0 at first station"); ax.set_ylabel("concentration C_i"); ax.set_ylim(0, 0.5)
    panel_label(ax, "a")

    # ternary
    ax = axs[1]; ax.set_aspect("equal"); ax.axis("off"); ax.grid(False)
    def tern(f):   # f: (n,3) fractions -> xy
        f = np.asarray(f); return np.c_[f[:, 1] + 0.5 * f[:, 2], (np.sqrt(3) / 2) * f[:, 2]]
    tri = tern(np.eye(3)); tri = np.vstack([tri, tri[0]])
    ax.plot(tri[:, 0], tri[:, 1], color=INK2, lw=0.8)
    for lab, (px, py), ha in zip(["A", "B", "C"], tern(np.eye(3)), ["right", "left", "center"]):
        ax.text(px + (-0.02 if ha == "right" else 0.02 if ha == "left" else 0), py + (0.02 if ha == "center" else 0), lab, ha=ha, va="bottom" if ha == "center" else "center", fontweight="bold")
    xs = np.linspace(-60, 0, 400)
    src = C0 / C0.sum(); sxy = tern(src[None, :])[0]
    acc = rej = 0
    for k in range(40):
        lam_k = lam * rng.uniform(0.6, 1.4, 3); c0_k = C0 * np.exp(-lam * L) * rng.uniform(0.7, 1.3, 3)   # draws around the "fitted" values at x=0
        Ck = c0_k[:, None] * np.exp(-lam_k[:, None] * xs[None, :]); f = (Ck / Ck.sum(0)).T
        d = np.sqrt(((f[:, :2] - src[:2]) ** 2).sum(1)); j = d.argmin()
        ok = d[j] < 0.05 and 0 < j < len(xs) - 1
        xy = tern(f)
        ax.plot(xy[:, 0], xy[:, 1], color=BLUE if ok else "#b9b8b2", lw=1.0 if ok else 0.6, alpha=0.9 if ok else 0.7, zorder=3 if ok else 2)
        acc += ok; rej += (not ok)
    ax.add_patch(plt.Circle(sxy, 0.05, fill=False, color=INK, lw=1.0, zorder=5))
    ax.plot(*sxy, "o", color=INK, ms=4, zorder=6); ax.text(sxy[0] + 0.07, sxy[1] - 0.02, "source\ncomposition\nand tolerance", fontsize=7.5, va="top")
    ax.text(0.5, -0.12, f"{acc} accepted (blue), {rej} rejected (grey) of 40 trials", transform=ax.transAxes, ha="center", fontsize=7.5, color=INK2)
    panel_label(ax, "b")

    ax = axs[2]
    # histogram from a larger draw
    xstar = []
    for k in range(2000):
        lam_k = lam * rng.uniform(0.6, 1.4, 3); c0_k = C0 * np.exp(-lam * L) * rng.uniform(0.7, 1.3, 3)
        Ck = c0_k[:, None] * np.exp(-lam_k[:, None] * xs[None, :]); f = (Ck / Ck.sum(0)).T
        d = np.sqrt(((f[:, :2] - src[:2]) ** 2).sum(1)); j = d.argmin()
        if d[j] < 0.05 and 0 < j < len(xs) - 1: xstar.append(-xs[j])
    xstar = np.array(xstar)
    ax.hist(xstar, bins=np.arange(0, 62, 2), color=BLUE, edgecolor="white", lw=0.5)
    ax.axvline(L, color=INK, lw=1.0); ax.text(L + 1, ax.get_ylim()[1] * 0.95, "true source", va="top", fontsize=7.5)
    ax.set_xlabel("distance upstream of first station (m)"); ax.set_ylabel("accepted trials")
    panel_label(ax, "c")
    fig.tight_layout(); save(fig, "fig1_concept")


# ============================================================== shared loaders
def load_fc():
    size = pd.read_csv(f"{res}/app_inputs/FortCuster_sampling_perSize.csv"); size["Filter"] = size["Filter"].astype(str)
    return size


SNYDER_VDEP = {  # Snyder 2024 Bayesian Vdep (m/d): season -> (total, 1.2, 0.4)
    "S23": (0.53, 0.64, 0.31), "Sm23": (0.73, 0.99, 0.47), "F22": (0.72, 0.91, 0.43), "W23": (0.33, 0.77, 0.37)}
SEASON_LAB = {"F22": "Fall 2022", "W23": "Winter 2023", "S23": "Spring 2023", "Sm23": "Summer 2023"}
SEASON_ORDER = ["F22", "W23", "S23", "Sm23"]


# ============================================================== Figure 3: Fort Custer profiles + Vdep comparison
def fig3():
    size = load_fc()
    d = O.normalise(size)
    reg = pd.read_csv(f"{res}/stage2_regressions_original_method.csv")
    reg = reg[reg.experiment == "FortCuster"]; reg["component"] = reg.component.astype(str)
    fig = plt.figure(figsize=(10.5, 5.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.35, 1], hspace=0.5, wspace=0.3)
    for j, se in enumerate(SEASON_ORDER):
        ax = fig.add_subplot(gs[0, j])
        g = d[d.Treatment == se]
        dn = O.denorm_factor(d, se)
        for comp in ["1.2", "0.4"]:
            gg = g[g.Filter == comp]
            xx = gg.Station_m if "Station_m" in gg else gg.Distance + 20
            ax.plot(gg.Distance + 20 + (1.5 if comp == "1.2" else -1.5), gg.Concentration.where(gg.Concentration > 0), "o", ms=3.5, color=COMP_COL[comp], alpha=0.75, mec="white", mew=0.4)
            r = reg[(reg.dataset == se) & (reg.component == comp)].iloc[0]
            xl = np.linspace(0, 80, 50)
            ax.plot(xl + 20, dn * np.exp(r.intercept + r.slope * xl), color=COMP_COL[comp], lw=1.6)
            # 95% band of the mean line from the pooled regression
            sub = gg[gg.ln_norm_c.notna()]; n = len(sub); xm = sub.Distance.mean(); sxx = ((sub.Distance - xm) ** 2).sum()
            resid = sub.ln_norm_c - (r.intercept + r.slope * sub.Distance); s2 = (resid ** 2).sum() / (n - 2)
            se_line = np.sqrt(s2 * (1 / n + (xl - xm) ** 2 / sxx)); t = stats.t.ppf(0.975, n - 2)
            ax.fill_between(xl + 20, dn * np.exp(r.intercept + r.slope * xl - t * se_line), dn * np.exp(r.intercept + r.slope * xl + t * se_line), color=COMP_COL[comp], alpha=0.12, lw=0)
        ax.set_yscale("log"); ax.set_ylim(1.5e-3, 6e-2); ax.set_yticks([2e-3, 5e-3, 1e-2, 2e-2, 5e-2]); ax.set_yticklabels(["0.002", "0.005", "0.01", "0.02", "0.05"]); ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); ax.set_title(SEASON_LAB[se]); ax.set_xticks([20, 40, 60, 80, 100])
        if j > 0: ax.set_yticklabels([])
        ax.set_xlabel("distance from injection (m)")
        if j == 0: ax.set_ylabel("eDNA (copies/mL per uS/cm)"); panel_label(ax, "a")
    handles = [Line2D([], [], color=COMP_COL[c], marker="o", ms=4, label=COMP_LAB[c]) for c in ["1.2", "0.4"]]
    fig.axes[0].legend(handles=handles, loc="lower left", frameon=False)

    ax = fig.add_subplot(gs[1, :2])
    app = reg.pivot(index="dataset", columns="component", values="Vdep_m_per_day")
    for comp, idx in [("total", 0), ("1.2", 1), ("0.4", 2)]:
        xs_ = [SNYDER_VDEP[s][idx] for s in SEASON_ORDER]; ys_ = [app.loc[s, comp] for s in SEASON_ORDER]
        ax.plot(xs_, ys_, "o", color=COMP_COL[comp], ms=6, mec="white", mew=0.6, label=COMP_LAB[comp])
        for s, xv, yv in zip(SEASON_ORDER, xs_, ys_):
            ax.annotate({"F22": "fall", "W23": "winter", "S23": "spring", "Sm23": "summer"}[s], (xv, yv), xytext=(5, -3), textcoords="offset points", fontsize=6.5, color=INK2)
    lim = (0, 1.15); ax.plot(lim, lim, color=INK2, lw=0.8, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("V_dep, Bayesian hierarchical fit, Snyder 2024 (m/d)"); ax.set_ylabel("V_dep, application regression (m/d)")
    ax.legend(frameon=False, loc="upper left"); ax.text(-0.16, 1.04, "b", transform=ax.transAxes, fontweight="bold", fontsize=10, va="bottom", ha="right")

    ax = fig.add_subplot(gs[1, 2:])
    w = 0.26; xs_ = np.arange(4)
    for k, comp in enumerate(["total", "1.2", "0.4"]):
        vals = [reg[(reg.dataset == s) & (reg.component == comp)]["lambda"].iloc[0] for s in SEASON_ORDER]
        lo = [-reg[(reg.dataset == s) & (reg.component == comp)]["slope_high"].iloc[0] for s in SEASON_ORDER]
        hi = [-reg[(reg.dataset == s) & (reg.component == comp)]["slope_low"].iloc[0] for s in SEASON_ORDER]
        ax.bar(xs_ + (k - 1) * w, vals, width=w - 0.03, color=COMP_COL[comp], label=COMP_LAB[comp])
        ax.errorbar(xs_ + (k - 1) * w, vals, yerr=[np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)], fmt="none", ecolor=INK, elinewidth=0.8, capsize=2)
    ax.axhline(0, color=INK2, lw=0.8); ax.set_xticks(xs_); ax.set_xticklabels([SEASON_LAB[s].split()[0] for s in SEASON_ORDER])
    ax.set_ylabel("removal rate lambda (1/m)"); ax.legend(frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.12)); panel_label(ax, "c")
    save(fig, "fig3_fortcuster_rates")


# ============================================================== Figure 4: SMIM r_h
PHASE_ORDER = ["S", "LB", "HB", "RS"]; PHASE_LAB = {"S": "early\nshaded", "LB": "early\nopen", "HB": "late\nopen", "RS": "late\nshaded"}
SUB_MARK = {"COB": "o", "PG": "s", "SAND": "^", "MIX": "D"}; SUB_LAB = {"COB": "cobble", "PG": "pea gravel", "SAND": "sand", "MIX": "mix"}


def fig4():
    S = pd.read_csv(f"{res}/stage3_smim_rh.csv"); S["component"] = S.component.astype(str)
    ok = S[S.k_CI_excludes_0 & S.rh_per_s.notna()].copy()
    ok["phase"] = ok.dataset.str.split("_").str[0]; ok["sub"] = ok.dataset.str.split("_").str[1]
    fig = plt.figure(figsize=(10.5, 4.2)); gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.2, 1], wspace=0.35)
    for j, (sp, lab) in enumerate([("CARP", "Common carp"), ("STLHD", "Steelhead trout")]):
        ax = fig.add_subplot(gs[0, j]); g = ok[ok.target == sp]
        for k, comp in enumerate(["1.2", "0.4", "total"]):
            gg = g[g.component == comp]
            for _, r in gg.iterrows():
                xpos = PHASE_ORDER.index(r.phase) + (k - 1) * 0.22 + np.random.default_rng(hash(r.dataset) % 1000).uniform(-0.05, 0.05)
                ax.errorbar(xpos, r.rh_per_s, yerr=[[r.rh_per_s - r.rh_at_k_CI_low], [r.rh_at_k_CI_high - r.rh_per_s]], fmt="none", ecolor=COMP_COL[comp], elinewidth=0.7, alpha=0.5)
                ax.plot(xpos, r.rh_per_s, SUB_MARK[r["sub"]], color=COMP_COL[comp], ms=5, mec="white", mew=0.5)
        ax.set_yscale("log"); ax.set_ylim(1e-5, 1e-1); ax.set_xticks(range(4)); ax.set_xticklabels([PHASE_LAB[p] for p in PHASE_ORDER])
        ax.set_title(f"{lab} ({len(g)} profiles)")
        rs_med = g.rs_per_s.median(); ax.axhline(rs_med, color=INK2, lw=0.8, ls="--"); ax.text(3.45, rs_med * 1.4, "r_s, water column (mesocosm fit)", fontsize=7, color=INK2, ha="right")
        if j == 0: ax.set_ylabel("hyporheic removal rate r_h (1/s)"); panel_label(ax, "a")
        else: panel_label(ax, "b")
    h1 = [Line2D([], [], color=COMP_COL[c], marker="o", ls="", ms=5, label=COMP_LAB[c]) for c in ["1.2", "0.4", "total"]]
    h2 = [Line2D([], [], color=INK2, marker=SUB_MARK[s], ls="", ms=5, label=SUB_LAB[s]) for s in SUB_MARK]
    fig.axes[0].legend(handles=h1 + h2, frameon=False, loc="upper right", ncol=2, fontsize=7)
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(ok.rh_per_s, ok.rh_t2_1e3, "o", color=BLUE, ms=4, mec="white", mew=0.4, label="t_2 = 1e3 s instead of 1e5 s")
    ax.plot(ok.rh_per_s, ok.rh_rs_0, "s", color=ORANGE, ms=4, mec="white", mew=0.4, label="r_s = 0 instead of fitted")
    lim = (1e-3, 5e-2); ax.plot(lim, lim, color=INK2, lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    for a_ in (ax.xaxis, ax.yaxis): a_.set_major_locator(matplotlib.ticker.FixedLocator([1e-3, 3e-3, 1e-2, 3e-2])); a_.set_major_formatter(matplotlib.ticker.FixedFormatter(["0.001", "0.003", "0.01", "0.03"])); a_.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("r_h, application settings (1/s)"); ax.set_ylabel("r_h, alternative setting (1/s)"); ax.legend(frameon=False, loc="upper left", fontsize=7); panel_label(ax, "c")
    save(fig, "fig4_smim_rh")


# ============================================================== reruns for histograms
def run_reference_station(size, spacing, treatment, nsim=1000, tmax=100, dt=2, tol=0.05, seed=1):
    g = size[size.Treatment == treatment]; smin = g.Station_m.min() if "Station_m" in g else None
    g = g.assign(Station_m=g.Distance + spacing)  # app inputs have Distance = station - first station
    src = g[g.Station_m == g.Station_m.min()]; src = src[src.Concentration > 0]
    nf = src.groupby("Replicate")["Filter"].transform("nunique"); src = src[nf == g.Filter.nunique()]
    samp = g[g.Station_m > g.Station_m.min()].copy(); samp["Distance"] = samp.Station_m - (g.Station_m.min() + spacing)
    d = O.normalise(samp); tab = O.table_var(d, treatment, by_size=True); f0 = O.source_composition(src)
    mc = O.mc_new(tab, nsim, dt, tmax, seed=seed)
    sol = O.df_close(mc, f0["mean"].to_numpy(), tol, tmax, denorm=O.denorm_factor(d, treatment))
    return sol, mc, f0


# ============================================================== Figure 5: source localisation, field
def fig5():
    size = load_fc()
    P = pd.read_csv(f"{res}/stage4a_reference_station_prediction.csv"); P = P[P.tol_rule == "app_0.05"]
    fig = plt.figure(figsize=(10.5, 6.2)); gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.25], hspace=0.55, wspace=0.35)
    mcs = {}
    for j, se in enumerate(SEASON_ORDER):
        ax = fig.add_subplot(gs[0, j]); sol, mc, f0 = run_reference_station(size, 20, se); mcs[se] = (mc, f0)
        ax.hist(sol.distance, bins=np.arange(0, 102, 2), color=BLUE, edgecolor="white", lw=0.4)
        ax.axvline(20, color=INK, lw=1.0)
        r = P[(P.experiment == "FortCuster") & (P.dataset == se)].iloc[0]
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18); y = ax.get_ylim()[1]
        ax.plot([r.x_P10, r.x_P90], [y * 0.92] * 2, color=ORANGE, lw=2); ax.plot(r.x_P50, y * 0.92, "o", color=ORANGE, ms=4)
        ax.set_title(f"{SEASON_LAB[se]}  (n = {len(sol)})"); ax.set_xlabel("x* upstream of 2nd station (m)")
        if j == 0: ax.set_ylabel("accepted trials"); panel_label(ax, "a")
    fig.axes[0].text(21, fig.axes[0].get_ylim()[1] * 0.6, "truth\n20 m", fontsize=7)
    fig.axes[0].text(2, fig.axes[0].get_ylim()[1] * 0.99, "P10 to P90, P50", color=ORANGE, fontsize=7, va="top")

    ax = fig.add_subplot(gs[1, :2]); mc, f0 = mcs["W23"]
    xgrid = -mc["time"]; f1 = mc["frac"][:, 0, :]
    for k in range(200):
        ax.plot(xgrid, f1[k], color="#b9b8b2", lw=0.5, alpha=0.6)
    for k in range(200):
        d = np.abs(f1[k] - f0["mean"].iloc[0]); j = d.argmin()
        if 0 < j < len(xgrid) - 1 and d[j] < 1e-9 or (0 < j < len(xgrid) - 1 and d[j] <= 0.05):
            ax.plot(xgrid[j], f1[k][j], "o", color=BLUE, ms=2.5, alpha=0.8)
    ax.axhline(f0["mean"].iloc[0], color=INK, lw=1.0); ax.text(99, f0["mean"].iloc[0] + 0.01, "first-station fraction f_0.4", ha="right", fontsize=7.5)
    ax.axhspan(f0["mean"].iloc[0] - 0.05, f0["mean"].iloc[0] + 0.05, color=ORANGE, alpha=0.12, lw=0); ax.text(99, f0["mean"].iloc[0] - 0.08, "tolerance band (no effect)", ha="right", fontsize=7.5, color=ORANGE)
    ax.set_xlabel("distance upstream of 2nd station (m)"); ax.set_ylabel("fraction of eDNA in 0.4 um class"); ax.set_title("Winter 2023: 200 trajectories and their crossings")
    ax.set_xlim(0, 100); panel_label(ax, "b")

    ax = fig.add_subplot(gs[1, 2:])
    Q = P[P.n_sol > 1].copy(); Q["truth"] = Q.known_x_m
    Q = Q.sort_values(["experiment", "x_P50"]).reset_index(drop=True)
    for i, r in Q.iterrows():
        col = BLUE if r.experiment == "LEEF2021" else ORANGE
        ax.plot([r.x_P10, r.x_P90], [i, i], color=col, lw=1.6, alpha=0.9 if r.covered else 0.45)
        ax.plot(r.x_P50, i, "o", color=col, ms=3)
        ax.plot(r.truth, i, "|", color=INK, ms=6, mew=1.2)
    ax.set_yticks([]); ax.set_xlim(0, 100); ax.set_xlabel("x* interval P10 to P90 and P50 (m); black tick = truth")
    ax.set_ylabel("34 field datasets, sorted by P50")
    ax.legend(handles=[Line2D([], [], color=BLUE, lw=2, label="ND-LEEF 2021 (truth 10 m)"), Line2D([], [], color=ORANGE, lw=2, label="Fort Custer (truth 20 m)"),
                       Line2D([], [], color=BLUE, lw=2, alpha=0.45, label="interval misses the truth")], frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7)
    panel_label(ax, "c")
    save(fig, "fig5_source_localisation_field")


# ============================================================== Figure 6: mesocosm validation recomputed + summary
def fig6():
    samp = pd.read_csv(f"{res}/app_inputs/Mesocosm_sampling_perSize_hours.csv"); samp["Filter"] = samp.Filter.astype(str)
    src = pd.read_csv(f"{res}/app_inputs/Mesocosm_source.csv"); src["Filter"] = src.Filter.astype(str)
    R = pd.read_csv(f"{res}/stage0b_mesocosm_table42_recomputed.csv"); A = R[R.run == "A_thesis_f_eps"].set_index("treatment")
    dn = O.normalise(samp)
    order = ["Dark_biof_Short", "Dark_biof_Long", "Dark_nosubs_Short", "Dark_nosubs_Long", "Dark_subs_Short", "Dark_subs_Long",
             "Light_biof_Short", "Light_biof_Long", "Light_nosubs_Short", "Light_nosubs_Long", "Light_subs_Short", "Light_subs_Long"]
    fig = plt.figure(figsize=(10.5, 9.0)); gs = fig.add_gridspec(4, 4, hspace=0.75, wspace=0.3, height_ratios=[1, 1, 1, 2.3])
    for i, tr in enumerate(order):
        ax = fig.add_subplot(gs[i // 4, i % 4]); a = A.loc[tr]
        tab = O.table_var(dn, tr, by_size=True); mc = O.mc_new(tab, 2000, 0.5, 100, seed=1)
        sol = O.df_close(mc, np.array([a.f1, a.f2, a.f3]), a.eps, 100, denorm=O.denorm_factor(dn, tr))
        ax.hist(sol.distance, bins=np.arange(0, 101, 2.5), color=BLUE, edgecolor="white", lw=0.3)
        ax.axvline(5.33, color=INK, lw=1.0); y = ax.get_ylim()[1]
        ax.plot([a.thesis_P10, a.thesis_P90], [y * 0.9] * 2, color=ORANGE, lw=2); ax.plot(a.thesis_P50, y * 0.9, "o", color=ORANGE, ms=3.5)
        ax.set_title(f"{a.acronym}  (n = {len(sol)})", fontsize=8.5); ax.set_xlim(0, 100)
        if i % 4 == 0: ax.set_ylabel("trials")
        if i >= 8: ax.set_xlabel("source distance (h)")
    fig.axes[0].text(7, fig.axes[0].get_ylim()[1] * 0.55, "truth 5.33 h", fontsize=7)
    fig.axes[0].text(50, fig.axes[0].get_ylim()[1] * 0.98, "thesis P10 to P90, P50", color=ORANGE, fontsize=7, va="top", ha="center")
    panel_label(fig.axes[0], "a")

    ax = fig.add_subplot(gs[3, :])
    P = pd.read_csv(f"{res}/stage4a_reference_station_prediction.csv"); P = P[(P.tol_rule == "app_0.05") & (P.n_sol > 1)]
    rows = [("3 components, mesocosm (h / 100 h)", A.x_P10 / 100, A.x_P50 / 100, A.x_P90 / 100, 5.33 / 100, AQUA)]
    L = P[P.experiment == "LEEF2021"]; F = P[P.experiment == "FortCuster"]
    rows += [("2 components, ND-LEEF (m / 100 m)", L.x_P10 / 100, L.x_P50 / 100, L.x_P90 / 100, 0.10, BLUE),
             ("2 components, Fort Custer (m / 100 m)", F.x_P10 / 100, F.x_P50 / 100, F.x_P90 / 100, 0.20, ORANGE)]
    y0 = 0; yt = []
    for lab, p10, p50, p90, truth, col in rows:
        idx = np.argsort(p50.to_numpy())
        for k, ii in enumerate(idx):
            yy = y0 + k
            ax.plot([p10.iloc[ii], p90.iloc[ii]], [yy, yy], color=col, lw=1.0, alpha=0.9); ax.plot(p50.iloc[ii], yy, "o", color=col, ms=2.2)
        ax.plot([truth, truth], [y0 - 0.5, y0 + len(idx) - 0.5], color=INK, lw=1.0)
        yt.append((y0 + len(idx) / 2, lab)); y0 += len(idx) + 2
    ax.set_yticks([t[0] for t in yt]); ax.set_yticklabels([t[1] for t in yt], fontsize=8); ax.set_xlim(0, 1)
    ax.set_xlabel("source-distance interval P10 to P90 and P50, as a fraction of the search window; black line = truth")
    panel_label(ax, "b")
    save(fig, "fig6_validation_vs_field")


if __name__ == "__main__":
    which = sys.argv[4] if len(sys.argv) > 4 else "all"
    for name, fn in [("fig1", fig1), ("fig3", fig3), ("fig4", fig4), ("fig5", fig5), ("fig6", fig6)]:
        if which in ("all", name): fn()
