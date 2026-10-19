"""
eDNA analysis tool, Shiny for Python edition.

Same four analyses as the R Shiny application (Hallack 2025, Ch. 4; github.com/dhallackla/edna_app),
driven by the Python port of its model routines in original.py, so every number the R app gives
can be checked here and vice versa. Input CSVs use the R app's templates:
  sampling: Treatment, Replicate, Filter, Target, Concentration, Distance  (or Hours)
  source:   Treatment, Replicate, Filter, Target, Concentration

Run:  python -m shiny run --port 7799 edna_py/app.py   (from the repository root)
"""
from __future__ import annotations
import io, os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shiny import App, reactive, render, ui

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import original as O

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(os.path.dirname(HERE), "app_R")
if not os.path.isdir(TEMPLATE_DIR):            # shinylive export: templates copied beside app.py
    TEMPLATE_DIR = os.path.join(HERE, "templates")
def sid(f):
    """input id safe for Shiny: letters, digits, underscore"""
    return re.sub(r"[^A-Za-z0-9_]", "_", str(f))

COLS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]

ABOUT = """
### eDNA analysis tool (Python edition)

Interprets eDNA concentration profiles from steady-state additions in streams, or time series in closed
systems. Data may be **total** or **by component** (particle size, amplicon length, any partition that decays
independently) and a function of **distance** or **time**. Four analyses, available depending on the data:

1. **LOD-distance**: furthest point at which eDNA is still above a limit of detection.
2. **SMIM-rates**: water-column vs hyporheic removal rates from the reactive stochastic mobile-immobile model
   at steady state (distance data only; needs U, D, Lambda, beta from a tracer test).
3. **Prediction**: distance upstream to a source of known composition and the concentration released there
   (by-component data only).
4. **Sensitivity**: composition upstream for chosen rates and distance (by-component data only).

Assumptions: independent first-order decay per component; a fixed source at steady state; conservative SMIM
parameters known. Conventions follow the R application exactly: replicate-maximum normalisation before
regression, non-detects dropped, slopes clipped at zero in the Monte Carlo, Euclidean distance over the first
n-1 fractions, solutions at the search boundaries discarded, histogram mode over 15 bins.

Models: Hallack (2025) PhD thesis, Univ. of Notre Dame, Ch. 4; Brandão-Dias et al. (2023) Mol. Ecol. Resour.
23:756; Badilla-Aguilar et al. (2024) Water Res. 259:121833. R application: github.com/dhallackla/edna_app.
"""

# ----------------------------------------------------------------------------- UI
def num(id_, label, value, **kw):
    return ui.input_numeric(id_, label, value, **kw)

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h5("Data"),
        ui.input_radio_buttons("data_type", "eDNA concentration", {"2": "by component", "1": "total"}, selected="2", inline=True),
        ui.input_radio_buttons("data_in", "Independent variable", {"1": "distance", "2": "time"}, selected="1", inline=True),
        ui.input_file("file1", "Sampling CSV", accept=[".csv"]),
        ui.output_ui("treat_ui"),
        ui.input_file("file2", "Source CSV (Prediction)", accept=[".csv"]),
        ui.output_ui("source_ui"),
        ui.hr(),
        ui.download_button("dl_samp", "Sampling template", class_="btn-sm"),
        ui.download_button("dl_src", "Source template", class_="btn-sm"),
        width=300,
    ),
    ui.navset_tab(
        ui.nav_panel("About", ui.markdown(ABOUT)),
        ui.nav_panel("Sampling Data",
            ui.layout_columns(
                ui.card(ui.card_header("Regressions (ln of normalised concentration vs x)"), ui.output_table("reg_table")),
                ui.card(ui.card_header("Data and fits"), ui.output_plot("reg_plot", height="420px")),
                col_widths=[5, 7])),
        ui.nav_panel("Source Data",
            ui.layout_columns(
                ui.card(ui.card_header("Source composition (mean and SD over replicates)"), ui.output_table("src_table")),
                ui.card(ui.card_header("Composition at the first station (mean and SD over replicates)"), ui.output_table("first_table")),
                col_widths=[6, 6])),
        ui.nav_panel("Prediction",
            ui.layout_columns(
                ui.card(ui.card_header("Settings"),
                    num("nsim", "Number of extrapolation trials", 1000, min=10, max=100000),
                    num("tmax", "Maximum distance (or time) upstream", 100, min=1),
                    num("dt", "Step", 2, min=0.01),
                    num("tol", "Tolerance epsilon (composition distance)", 0.05, min=0.001, max=0.99, step=0.01),
                    num("d_max", "Max admissible distance (blank = none)", None),
                    num("c0_max", "Max admissible released C0 (blank = none)", None),
                    num("seed", "Random seed", 1),
                    ui.input_action_button("submit", "Run", class_="btn-primary")),
                ui.card(ui.card_header("Result"), ui.output_code("pred_text"), ui.output_plot("pred_hist", height="300px")),
                col_widths=[3, 9]),
            ui.layout_columns(
                ui.card(ui.card_header("Back-extrapolations (accepted in colour)"), ui.output_plot("pred_back", height="360px")),
                ui.card(ui.card_header("Composition trajectories and source"), ui.output_plot("pred_traj", height="360px")),
                col_widths=[6, 6]),
            ui.card(ui.card_header("Accepted trials (first 15)"), ui.output_table("pred_table"))),
        ui.nav_panel("Sensitivity",
            ui.layout_columns(
                ui.card(ui.card_header("Rates and intercepts (initialised to the regression, bounded by its 95% CI)"), ui.output_ui("sens_ui")),
                ui.card(ui.card_header("Composition upstream"), ui.output_plot("sens_plot", height="420px")),
                col_widths=[4, 8])),
        ui.nav_panel("LOD",
            ui.layout_columns(
                ui.card(ui.card_header("Limit of detection"), num("LOD", "LOD (copies/mL)", 0.05, min=1e-6), ui.output_table("lod_table")),
                ui.card(ui.card_header("Fitted profiles and LOD"), ui.output_plot("lod_plot", height="420px")),
                col_widths=[4, 8])),
        ui.nav_panel("SMIM",
            ui.layout_columns(
                ui.card(ui.card_header("Conservative parameters (tracer test) and water-column rate"),
                    ui.output_ui("smim_filter_ui"),
                    num("Us", "U (m/s)", 0.0982), num("Ds", "D (m^2/s)", 0.0345), num("lambdas", "Lambda (1/s)", 0.085),
                    num("betas", "beta", 0.898), num("logT1", "log10 t1 (s)  [R app: 0]", 0), num("logT2", "log10 t2 (s)  [R app: 5]", 5),
                    num("rss", "r_s, water column (1/s)", 3e-5)),
                ui.card(ui.card_header("Inversion for r_h"), ui.output_code("smim_text"),
                    ui.card_header("Forward profile for chosen rates"),
                    ui.layout_columns(num("rh_fwd", "r_h (1/s)", 1e-3), num("rs_fwd", "r_s (1/s)", 3e-5), num("xlim", "reach length (m)", 100), num("xloc", "report C at x (m)", 50), col_widths=[3, 3, 3, 3]),
                    ui.output_plot("smim_plot", height="320px")),
                col_widths=[4, 8])),
        id="tabs",
    ),
    title="eDNA analysis tool (Python edition)",
)


# ----------------------------------------------------------------------------- server
def server(input, output, session):

    @reactive.calc
    def raw():
        f = input.file1()
        if not f:
            return None
        d = pd.read_csv(f[0]["datapath"], encoding="utf-8-sig")
        if input.data_in() == "2" and "Hours" in d.columns:
            d = d.rename(columns={"Hours": "Distance"})
        if "Filter" in d.columns:
            d["Filter"] = d["Filter"].astype(str)
        return d

    @reactive.calc
    def by_size():
        return input.data_type() == "2"

    @output
    @render.ui
    def treat_ui():
        d = raw()
        if d is None:
            return ui.help_text("Upload a sampling CSV to begin.")
        return ui.TagList(
            ui.input_select("treat_var", "Treatment", choices=sorted(d["Treatment"].astype(str).unique()) + ["ALL"]),
            ui.input_select("target_var", "Target", choices=sorted(d["Target"].astype(str).unique())) if "Target" in d else None,
        )

    @reactive.calc
    def data_n():
        d = raw()
        if d is None or not input.treat_var():
            return None
        if "Target" in d.columns and input.target_var():
            d = d[d["Target"].astype(str) == input.target_var()]
        if input.treat_var() == "ALL":
            d = d.assign(Treatment="ALL")
        if not by_size():
            d = d.drop(columns=[c for c in ["Filter"] if c in d.columns])
        return O.normalise(d)

    @reactive.calc
    def treatment():
        return input.treat_var() if input.treat_var() else None

    @reactive.calc
    def table():
        d = data_n()
        if d is None:
            return None
        return O.table_var(d, treatment(), by_size=by_size())

    @reactive.calc
    def filters():
        t = table()
        return [] if t is None else list(t["Filter"].unique())

    # ---- Sampling Data
    @output
    @render.table
    def reg_table():
        t = table()
        if t is None:
            return pd.DataFrame({"note": ["upload data and choose a treatment"]})
        out = t[["Filter", "term", "estimate", "conf.low", "conf.high", "n", "r.squared"]].copy()
        out["term"] = out["term"].replace({"Distance": "slope (-lambda)", "(Intercept)": "intercept"})
        return out.round(5)

    @output
    @render.plot
    def reg_plot():
        d = data_n(); t = table()
        fig, ax = plt.subplots()
        if d is None or t is None:
            ax.text(0.5, 0.5, "no data", ha="center"); return fig
        g = d[d["Treatment"] == treatment()]
        for i, f in enumerate(filters()):
            gg = g[g["Filter"] == f]
            ax.plot(gg["Distance"], gg["ln_norm_c"], "o", color=COLS[i % 6], ms=4, alpha=0.7, label=f"{f}")
            r = t[t["Filter"] == f]
            ic = r[r.term == "(Intercept)"]["estimate"].iloc[0]; sl = r[r.term == "Distance"]["estimate"].iloc[0]
            xx = np.linspace(g["Distance"].min(), g["Distance"].max(), 50)
            ax.plot(xx, ic + sl * xx, color=COLS[i % 6])
        ax.set_xlabel("distance (m)" if input.data_in() == "1" else "time (h)"); ax.set_ylabel("ln(normalised concentration)")
        ax.legend(title="component", frameon=False); ax.grid(alpha=0.3)
        return fig

    # ---- Source Data
    @reactive.calc
    def source_df():
        f = input.file2()
        if not f:
            return None
        s = pd.read_csv(f[0]["datapath"], encoding="utf-8-sig"); s["Filter"] = s["Filter"].astype(str)
        return s

    @output
    @render.ui
    def source_ui():
        s = source_df()
        if s is None:
            n = max(len(filters()), 2)
            return ui.TagList(ui.help_text("No source CSV: enter the source fractions by hand (first n-1 components)."),
                              *[num(f"fr0{i}", f"f_{filters()[i] if i < len(filters()) else i+1}", round(1 / n, 3), min=0, max=1, step=0.01) for i in range(n - 1)])
        choices = {"target": "by Target", "treat": "by Treatment"}
        return ui.TagList(
            ui.input_select("src_mode", "Select source rows", choices),
            ui.input_select("src_target", "Target", choices=sorted(s["Target"].astype(str).unique())) if "Target" in s else None,
            ui.input_select("src_treat", "Treatment", choices=sorted(s["Treatment"].astype(str).unique())),
        )

    @reactive.calc
    def f_source():
        """Filter, mean, std in the order of the sampling-data filters."""
        s = source_df(); fl = filters()
        if s is not None:
            try:
                sub = s[s["Target"].astype(str) == input.src_target()] if input.src_mode() == "target" and "Target" in s else s[s["Treatment"].astype(str) == input.src_treat()]
            except Exception:
                sub = s
            comp = O.source_composition(sub)
            comp = comp.set_index("Filter").reindex(fl).reset_index()
            return comp
        vals = []
        for i in range(len(fl) - 1):
            try:
                vals.append(float(input[f"fr0{i}"]()))
            except Exception:
                vals.append(1 / len(fl))
        vals.append(max(0.0, 1 - sum(vals)))
        return pd.DataFrame({"Filter": fl, "mean": vals, "std": np.nan})

    @output
    @render.table
    def src_table():
        if not by_size() or not filters():
            return pd.DataFrame({"note": ["needs by-component sampling data"]})
        return f_source().round(4)

    @output
    @render.table
    def first_table():
        d = data_n()
        if d is None or not by_size():
            return pd.DataFrame({"note": ["needs by-component sampling data"]})
        return O.first_sample_composition(d, treatment()).round(4)

    # ---- Prediction
    @reactive.calc
    @reactive.event(input.submit)
    def mc():
        t = table()
        if t is None or not by_size() or len(filters()) < 2:
            return None
        return O.mc_new(t, int(input.nsim()), float(input.dt()), float(input.tmax()), seed=int(input.seed()))

    @reactive.calc
    def solutions():
        m = mc()
        if m is None:
            return None
        f0 = f_source()["mean"].to_numpy(float)
        dn = O.denorm_factor(data_n(), treatment())
        return O.df_close(m, f0, float(input.tol()), float(input.tmax()), denorm=dn,
                          d_max=input.d_max() if input.d_max() else None, c0_max=input.c0_max() if input.c0_max() else None)

    @output
    @render.code
    def pred_text():
        sol = solutions()
        if sol is None:
            return "Load by-component data with at least two components, set the source composition, press Run."
        s = O.summarise_solutions(sol); unit = "m" if input.data_in() == "1" else "h"
        return (f"Accepted trials: {s['n_sol']} of {int(input.nsim())}\n"
                f"Source distance ({unit}):  mean {s['x_mean']:.2f}   SD {s['x_sd']:.2f}   mode {s['x_mode']:.2f}   P10/P50/P90 {s['x_P10']:.1f} / {s['x_P50']:.1f} / {s['x_P90']:.1f}\n"
                f"Source concentration log10(copies/mL):  mean {s['logC0_mean']:.2f}   SD {s['logC0_sd']:.2f}   mode {s['logC0_mode']:.2f}   P10/P50/P90 {s['logC0_P10']:.2f} / {s['logC0_P50']:.2f} / {s['logC0_P90']:.2f}")

    @output
    @render.plot
    def pred_hist():
        sol = solutions(); fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
        if sol is None or len(sol) == 0:
            axs[0].text(0.5, 0.5, "no solutions", ha="center"); return fig
        axs[0].hist(sol["distance"], bins=15, color=COLS[0], edgecolor="white"); axs[0].set_xlabel("source distance"); axs[0].set_ylabel("trials")
        axs[1].hist(np.log10(sol["C0denorm"]), bins=15, color=COLS[1], edgecolor="white"); axs[1].set_xlabel("log10 source concentration (copies/mL)")
        for a in axs: a.grid(alpha=0.3)
        fig.tight_layout(); return fig

    @output
    @render.plot
    def pred_back():
        m = mc(); sol = solutions(); d = data_n()
        fig, axs = plt.subplots(1, max(len(filters()), 1), figsize=(9, 3.4), sharey=True)
        axs = np.atleast_1d(axs)
        if m is None:
            axs[0].text(0.5, 0.5, "run the prediction", ha="center"); return fig
        acc = set(sol["Nsim"].to_numpy()) if sol is not None else set()
        g = d[d["Treatment"] == treatment()]
        for i, f in enumerate(m["filters"]):
            ax = axs[i]
            for k in range(min(m["ln_norm_c"].shape[0], 300)):
                ax.plot(m["time"], m["ln_norm_c"][k, i], color=(COLS[0] if (k + 1) in acc else "#c9c8c2"), lw=0.5, alpha=0.6, zorder=2 if (k + 1) in acc else 1)
            gg = g[g["Filter"] == f]; ax.plot(gg["Distance"], gg["ln_norm_c"], "o", color="k", ms=3.5, zorder=3)
            ax.set_title(f"component {f}"); ax.set_xlabel("x (negative = upstream)"); ax.grid(alpha=0.3)
        axs[0].set_ylabel("ln(normalised concentration)"); fig.tight_layout(); return fig

    @output
    @render.plot
    def pred_traj():
        m = mc(); sol = solutions(); f0 = f_source()["mean"].to_numpy(float) if filters() else None
        fig, ax = plt.subplots(figsize=(5, 4.2))
        if m is None:
            ax.text(0.5, 0.5, "run the prediction", ha="center"); return fig
        acc = set(sol["Nsim"].to_numpy()) if sol is not None else set()
        fr = m["frac"]; n = fr.shape[1]
        if n == 3:
            def tern(f): return np.c_[f[..., 1] + 0.5 * f[..., 2], (np.sqrt(3) / 2) * f[..., 2]]
            tri = tern(np.eye(3)); ax.plot(np.r_[tri[:, 0], tri[0, 0]], np.r_[tri[:, 1], tri[0, 1]], color="k", lw=0.8)
            for k in range(min(fr.shape[0], 300)):
                xy = tern(fr[k].T); ax.plot(xy[:, 0], xy[:, 1], color=(COLS[0] if (k + 1) in acc else "#c9c8c2"), lw=0.6, alpha=0.6)
            s = tern(f0 / f0.sum()); ax.plot(*s, "o", color="k"); ax.add_patch(plt.Circle(s, float(input.tol()), fill=False, color="k"))
            for lab, p in zip(m["filters"], tern(np.eye(3))): ax.text(p[0], p[1], lab, ha="center", va="bottom")
            ax.set_aspect("equal"); ax.axis("off")
        else:
            for k in range(min(fr.shape[0], 300)):
                ax.plot(-m["time"], fr[k, 0], color=(COLS[0] if (k + 1) in acc else "#c9c8c2"), lw=0.6, alpha=0.6)
            ax.axhline(f0[0], color="k"); ax.axhspan(f0[0] - float(input.tol()), f0[0] + float(input.tol()), color=COLS[1], alpha=0.15, lw=0)
            ax.set_xlabel("distance upstream"); ax.set_ylabel(f"fraction of component {m['filters'][0]}"); ax.grid(alpha=0.3)
            if n > 3: ax.set_title(f"{n} components: first fraction shown; distance uses the first {n-1}")
        return fig

    @output
    @render.table
    def pred_table():
        sol = solutions()
        if sol is None:
            return pd.DataFrame()
        return sol.head(15).round(4)

    # ---- Sensitivity
    @output
    @render.ui
    def sens_ui():
        t = table()
        if t is None or not by_size():
            return ui.help_text("needs by-component sampling data")
        items = []
        for f in filters():
            r = t[t["Filter"] == f]; sl = r[r.term == "Distance"].iloc[0]; ic = r[r.term == "(Intercept)"].iloc[0]
            lam, lo, hi = -sl["estimate"], max(-sl["conf.high"], 0), max(-sl["conf.low"], 1e-6)
            items.append(ui.input_slider(f"rem_{sid(f)}", f"lambda {f}", min=round(lo, 4), max=round(hi, 4), value=round(max(lam, 0), 4), step=round((hi - lo) / 100, 5) or 1e-4))
            items.append(ui.input_slider(f"int_{sid(f)}", f"intercept {f}", min=round(ic["conf.low"], 3), max=round(ic["conf.high"], 3), value=round(ic["estimate"], 3), step=0.01))
        items.append(ui.input_slider("time", "distance (or time) upstream of the first sample", min=0, max=100, value=10))
        return ui.TagList(*items)

    @output
    @render.plot
    def sens_plot():
        t = table(); fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
        if t is None or not by_size():
            axs[0].text(0.5, 0.5, "needs by-component data", ha="center"); return fig
        fl = filters()
        try:
            lam = np.array([float(input[f"rem_{sid(f)}"]()) for f in fl]); ic = np.array([float(input[f"int_{sid(f)}"]()) for f in fl]); X = float(input.time())
        except Exception:
            return fig
        xs = np.linspace(0, max(X, 1), 200)
        C = np.exp(ic[:, None] + lam[:, None] * xs[None, :]); fr = C / C.sum(0)
        for i, f in enumerate(fl):
            axs[0].plot(-xs, np.log(C[i]), color=COLS[i % 6], label=f); axs[1].plot(xs, fr[i], color=COLS[i % 6], label=f)
        axs[0].set_xlabel("x (negative = upstream)"); axs[0].set_ylabel("ln(normalised concentration)"); axs[0].legend(frameon=False)
        axs[1].set_xlabel("distance upstream"); axs[1].set_ylabel("fraction"); axs[1].set_ylim(0, 1)
        axs[1].set_title(f"composition at {X:g}: " + ", ".join(f"{v:.2f}" for v in fr[:, -1]))
        for a in axs: a.grid(alpha=0.3)
        fig.tight_layout(); return fig

    # ---- LOD
    @output
    @render.table
    def lod_table():
        t = table(); d = data_n()
        if t is None:
            return pd.DataFrame({"note": ["upload data"]})
        return O.lod_distance(t, O.max_factor(d, treatment()), float(input.LOD())).round(2)

    @output
    @render.plot
    def lod_plot():
        t = table(); d = data_n(); fig, ax = plt.subplots()
        if t is None:
            ax.text(0.5, 0.5, "upload data", ha="center"); return fig
        mf = O.max_factor(d, treatment()); lnmax = mf["lnmax"].mean(); L = O.lod_distance(t, mf, float(input.LOD()))
        xmax = np.nanmax(L[["avg", "upper", "lower"]].to_numpy(float)) if np.isfinite(L[["avg", "upper", "lower"]].to_numpy(float)).any() else d["Distance"].max() * 2
        xx = np.linspace(0, max(xmax * 1.05, d["Distance"].max()), 100)
        g = d[d["Treatment"] == treatment()]
        for i, f in enumerate(filters()):
            r = t[t["Filter"] == f]; ic = r[r.term == "(Intercept)"]["estimate"].iloc[0]; sl = r[r.term == "Distance"]["estimate"].iloc[0]
            ax.plot(xx, ic + sl * xx + lnmax, color=COLS[i % 6], label=f)
            gg = g[g["Filter"] == f]; ax.plot(gg["Distance"], gg["ln_c"], "o", color=COLS[i % 6], ms=3.5, alpha=0.7)
        ax.axhline(np.log(float(input.LOD())), color="red", ls="--", label="LOD")
        ax.set_xlabel("distance (m)" if input.data_in() == "1" else "time (h)"); ax.set_ylabel("ln(concentration, copies/mL)"); ax.legend(frameon=False); ax.grid(alpha=0.3)
        return fig

    # ---- SMIM
    @output
    @render.ui
    def smim_filter_ui():
        if by_size() and filters():
            return ui.input_select("size2smim", "Component for k", choices=filters())
        return ui.help_text("k taken from the total regression")

    @reactive.calc
    def k_obs():
        t = table()
        if t is None or input.data_in() == "2":
            return None
        f = input.size2smim() if (by_size() and "size2smim" in input) else (filters()[0] if filters() else None)
        if f is None:
            return None
        r = t[(t["Filter"] == f) & (t.term == "Distance")]
        return None if r.empty else -float(r["estimate"].iloc[0])

    @output
    @render.code
    def smim_text():
        k = k_obs()
        if k is None:
            return "Needs distance data (SMIM is not available for time series)."
        rh, res = O.fit_rh(k, float(input.rss()), float(input.Us()), float(input.Ds()), float(input.lambdas()), float(input.betas()), logT1=float(input.logT1()), logT2=float(input.logT2()))
        k0 = O.k_func(float(input.rss()), 0.0, float(input.Us()), float(input.Ds()), float(input.lambdas()), float(input.betas()), float(input.logT1()), float(input.logT2()))
        return (f"k from regression: {k:.5g} 1/m\n"
                f"r_h by least squares: {rh:.3e} 1/s   (|k_model - k| = {res:.2e})\n"
                f"k the model gives with r_h = 0 (water column only): {k0:.3e} 1/m\n"
                f"r_h / r_s = {rh / float(input.rss()):.1f}" if input.rss() else "")

    @output
    @render.plot
    def smim_plot():
        fig, ax = plt.subplots(figsize=(8, 3.4))
        U, D, L, B, lt1, lt2 = (float(input.Us()), float(input.Ds()), float(input.lambdas()), float(input.betas()), float(input.logT1()), float(input.logT2()))
        xs = np.linspace(0, float(input.xlim()), 200)
        k_fwd = O.k_func(float(input.rs_fwd()), float(input.rh_fwd()), U, D, L, B, lt1, lt2)
        ax.plot(xs, np.exp(-k_fwd * xs), color=COLS[1], label=f"chosen rates: k = {k_fwd:.4g} 1/m")
        k = k_obs()
        if k is not None:
            ax.plot(xs, np.exp(-k * xs), color=COLS[0], label=f"fitted: k = {k:.4g} 1/m")
        xl = float(input.xloc()); ax.axvline(xl, color="k", lw=0.8, ls=":"); ax.text(xl, 0.95, f"C/C0 at {xl:g} m = {np.exp(-k_fwd * xl):.3g}", fontsize=8)
        ax.set_xlabel("distance (m)"); ax.set_ylabel("C / C0"); ax.legend(frameon=False); ax.grid(alpha=0.3)
        return fig

    # ---- templates
    @render.download(filename="template_sampling_perSize.csv")
    def dl_samp():
        p = os.path.join(TEMPLATE_DIR, "template_sampling_perSize.csv")
        yield open(p, "rb").read() if os.path.exists(p) else b"Treatment,Replicate,Filter,Target,Concentration,Distance\n"

    @render.download(filename="template_source.csv")
    def dl_src():
        p = os.path.join(TEMPLATE_DIR, "template_source.csv")
        yield open(p, "rb").read() if os.path.exists(p) else b"Treatment,Replicate,Filter,Target,Concentration\n"


app = App(app_ui, server)
