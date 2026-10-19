"""
Python port of the model code in Daniel Hallack's original Shiny app
(github.com/dhallackla/edna_app, appv6+.R and aux_fun.R, commit 4ac13a7).

Every function names the R code it reproduces. Where the R code does something
the thesis text does not say, that is reproduced too (replicate-max normalisation,
n-1 fractions in the composition distance, boundary solutions discarded, histogram
mode, t1 = 1 s and t2 = 1e5 s). Two things cannot be reproduced bit-for-bit:
R's random stream (set.seed(1) + runif) and the SANN optimiser; both are replaced
by their deterministic equivalents and flagged in the docstrings.

Data frame conventions follow the app's CSV templates:
  sampling: Treatment, Replicate, Filter, Target, Concentration, Distance (or Hours)
  source:   Treatment, Replicate, Filter, Target, Concentration
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------
# Data preparation: raw.data.n() in appv6+.R
# --------------------------------------------------------------------------

def normalise(df: pd.DataFrame, x_col: str = "Distance") -> pd.DataFrame:
    """raw.data.n(): keep Distance >= 0; within each Treatment x Replicate divide
    every Concentration (all filters, all distances) by that replicate's maximum;
    ln of zero is NA. `frac` is the composition at each (Treatment, Replicate, x)."""
    d = df.copy()
    if x_col != "Distance":
        d = d.rename(columns={x_col: "Distance"})
    d["Distance"] = d["Distance"].astype(float)
    d["Concentration"] = pd.to_numeric(d["Concentration"], errors="coerce")
    d = d[d["Distance"] >= 0].copy()
    if "Filter" not in d.columns:          # total data: the app has no Filter column
        d["Filter"] = "total"
    d["Filter"] = d["Filter"].astype(str)
    g = d.groupby(["Treatment", "Replicate"])["Concentration"]
    d["max_c"] = g.transform("max")
    d["norm_c"] = d["Concentration"] / d["max_c"]
    d["ln_c"] = np.log(d["Concentration"].where(d["Concentration"] > 0))      # 0 or NA -> NaN
    d["ln_norm_c"] = np.log(d["norm_c"].where(d["norm_c"] > 0))
    s = d.groupby(["Treatment", "Replicate", "Distance"])["norm_c"].transform("sum")
    d["frac"] = d["norm_c"] / s
    return d


def filters_sorted(d: pd.DataFrame) -> list[str]:
    """dat()$sizes: filters sorted numerically (the regression table is sorted as a
    factor, i.e. as strings; the two orders agree for the filter names in use)."""
    return [str(f) for f in sorted(d["Filter"].unique(), key=lambda s: float(s))]


# --------------------------------------------------------------------------
# Regressions: table_var() in appv6+.R  (lm(ln_norm_c ~ Distance) + broom::tidy)
# --------------------------------------------------------------------------

def _lm(x, y, conf=0.95):
    """lm(y ~ x) with broom::tidy(conf.int = TRUE): t-based CIs on both terms."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]
    n = len(x)
    if n < 3 or len(np.unique(x)) < 2:
        return None
    r = stats.linregress(x, y)
    t = stats.t.ppf(0.5 + conf / 2, n - 2)
    return pd.DataFrame({
        "term": ["(Intercept)", "Distance"],
        "estimate": [r.intercept, r.slope],
        "std.error": [r.intercept_stderr, r.stderr],
        "p.value": [2 * stats.t.sf(abs(r.intercept / r.intercept_stderr), n - 2), r.pvalue],
        "conf.low": [r.intercept - t * r.intercept_stderr, r.slope - t * r.stderr],
        "conf.high": [r.intercept + t * r.intercept_stderr, r.slope + t * r.stderr],
        "n": [n, n],
        "r.squared": [r.rvalue ** 2] * 2,
    })


def table_var(d: pd.DataFrame, treatment: str | None = None, by_size: bool = True) -> pd.DataFrame:
    """One regression per Filter (by_size) or one for the treatment (total), on
    ln_norm_c vs Distance pooling replicates. treatment=None pools all ('ALL')."""
    r = d if treatment is None else d[d["Treatment"] == treatment]
    rows = []
    if by_size:
        for f in filters_sorted(r):
            sub = r[r["Filter"] == f]
            t = _lm(sub["Distance"], sub["ln_norm_c"])
            if t is not None:
                t.insert(0, "Filter", f); rows.append(t)
    else:
        t = _lm(r["Distance"], r["ln_norm_c"])
        if t is not None:
            t.insert(0, "Filter", "total"); rows.append(t)
    if not rows:
        return pd.DataFrame(columns=["Filter", "term", "estimate", "std.error", "p.value", "conf.low", "conf.high", "n", "r.squared"])
    return pd.concat(rows, ignore_index=True)


def denorm_factor(d: pd.DataFrame, treatment: str) -> float:
    """denorm_factor(): mean of the replicate maxima for the treatment."""
    return d.loc[d["Treatment"] == treatment, "max_c"].drop_duplicates().mean()


def max_factor(d: pd.DataFrame, treatment: str) -> pd.DataFrame:
    """max_factor_totalc(): one row per replicate with maxC and lnmax."""
    m = d[d["Treatment"] == treatment].groupby("Replicate")["max_c"].first().reset_index()
    m["lnmax"] = np.log(m["max_c"])
    return m.rename(columns={"max_c": "maxC"})


# --------------------------------------------------------------------------
# Source composition: fr0_species() in appv6+.R
# --------------------------------------------------------------------------

def source_composition(src: pd.DataFrame) -> pd.DataFrame:
    """Fraction per Treatment x Replicate, then mean and SD per Filter across
    replicates. Returns Filter, mean, std sorted by Filter."""
    s = src.copy()
    s["Filter"] = s["Filter"].astype(str)
    s["Concentration"] = pd.to_numeric(s["Concentration"], errors="coerce")
    s["frac"] = s["Concentration"] / s.groupby(["Treatment", "Replicate"])["Concentration"].transform("sum")
    out = s.groupby("Filter")["frac"].agg(mean="mean", std="std").reset_index()
    out = out.iloc[np.argsort([float(f) for f in out["Filter"]])].reset_index(drop=True)
    return out


def first_sample_composition(d: pd.DataFrame, treatment: str) -> pd.DataFrame:
    """f_1st_sample(): composition at Distance == 0 (frac < 1), mean and SD over
    replicates per Filter. This is what the app plots as the first-sample point."""
    r = d[(d["Treatment"] == treatment) & (d["Distance"] == 0) & (d["frac"] < 1)]
    out = r.groupby("Filter")["frac"].agg(mean="mean", std="std").reset_index()
    out = out.iloc[np.argsort([float(f) for f in out["Filter"]])].reset_index(drop=True)
    return out


# --------------------------------------------------------------------------
# Monte Carlo back-extrapolation: MC_new() + single_pred() in aux_fun.R
# --------------------------------------------------------------------------

def mc_new(table: pd.DataFrame, nsim: int, dt: float, tmax: float, seed: int | None = 1):
    """For each trial draw (intercept, slope) per Filter uniformly within their
    95% CIs, clip slope at 0 (min(slope, 0)), back-extrapolate over
    time = -seq(0, tmax, dt). Returns dict with time (nt,), filters, ln_norm_c
    (nsim, nf, nt), norm_c, frac, and the drawn parameters.
    R uses set.seed(1) + runif; numpy's generator differs, so trials are not the
    same draws, only the same distribution."""
    rng = np.random.default_rng(seed)
    filt = list(table["Filter"].unique())
    nf = len(filt)
    time = -np.arange(0, tmax + 1e-9, dt)
    nt = len(time)
    icpt = np.empty((nsim, nf)); slope = np.empty((nsim, nf))
    for i, f in enumerate(filt):
        ti = table[(table["Filter"] == f) & (table["term"] == "(Intercept)")].iloc[0]
        ts = table[(table["Filter"] == f) & (table["term"] == "Distance")].iloc[0]
        icpt[:, i] = rng.uniform(ti["conf.low"], ti["conf.high"], nsim)
        slope[:, i] = rng.uniform(ts["conf.low"], ts["conf.high"], nsim)
    slope_used = np.minimum(slope, 0.0)
    ln = icpt[:, :, None] + slope_used[:, :, None] * time[None, None, :]
    norm_c = np.exp(ln)
    csum = norm_c.sum(axis=1, keepdims=True)
    frac = norm_c / csum
    return dict(time=time, filters=filt, ln_norm_c=ln, norm_c=norm_c, frac=frac,
                intercept=icpt, slope=slope, slope_used=slope_used)


def df_close(mc: dict, f_source_mean, tol: float, tmax: float, denorm: float = 1.0,
             d_max: float | None = None, c0_max: float | None = None) -> pd.DataFrame:
    """df_close() + df_close_C0(): per trial, the grid point of minimum Euclidean
    distance between the trajectory's first n-1 fractions and the source's first
    n-1 fractions; accepted if distance <= tol; solutions at time == 0 and
    time == -tmax are discarded; optional d_max / c0_max cutoffs.
    Returns one row per accepted trial: Nsim, distance (= -time, positive
    upstream), SSe, C0 (normalised) and C0denorm."""
    frac = mc["frac"]; time = mc["time"]
    nf = frac.shape[1]
    f0 = np.asarray(f_source_mean, float)[: nf - 1]
    diff = frac[:, : nf - 1, :] - f0[None, :, None]
    sse = np.sqrt((diff ** 2).sum(axis=1))            # (nsim, nt)
    k = sse.argmin(axis=1)
    minsse = sse[np.arange(len(k)), k]
    tstar = time[k]
    ok = (minsse <= tol) & (tstar != 0) & (tstar != -tmax)
    c0 = mc["norm_c"][np.arange(len(k)), :, k].sum(axis=1)
    out = pd.DataFrame({"Nsim": np.arange(1, len(k) + 1), "distance": -tstar, "SSe": minsse,
                        "C0": c0, "C0denorm": c0 * denorm})[ok]
    if d_max is not None and d_max > 0:
        out = out[out["distance"] <= d_max]
    if c0_max is not None and c0_max > 0:
        out = out[out["C0denorm"] <= c0_max]
    return out.reset_index(drop=True)


def hist_mode(x, bins: int = 15) -> float:
    """mode() in aux_fun.R: centre of the fullest of 15 equal-width bins."""
    x = np.asarray(x, float)
    if len(x) <= 1 or x.max() == x.min():
        return float("nan")
    counts, edges = np.histogram(x, bins=np.linspace(x.min(), x.max(), bins + 1))
    i = counts.argmax()
    return float(0.5 * (edges[i] + edges[i + 1]))


def summarise_solutions(sol: pd.DataFrame) -> dict:
    """What the app prints (mean, SD, histogram mode) plus the P10/P50/P90 the
    thesis tables report (computed outside the app)."""
    x = sol["distance"].to_numpy(float)
    c = np.log10(sol["C0denorm"].to_numpy(float)) if len(sol) else np.array([])
    q = lambda v, p: float(np.quantile(v, p)) if len(v) else float("nan")
    return dict(n_sol=int(len(sol)),
                x_mean=float(x.mean()) if len(x) else np.nan, x_sd=float(x.std(ddof=1)) if len(x) > 1 else np.nan,
                x_mode=hist_mode(x), x_P10=q(x, .1), x_P50=q(x, .5), x_P90=q(x, .9),
                logC0_mean=float(c.mean()) if len(c) else np.nan, logC0_sd=float(c.std(ddof=1)) if len(c) > 1 else np.nan,
                logC0_mode=hist_mode(c), logC0_P10=q(c, .1), logC0_P50=q(c, .5), logC0_P90=q(c, .9))


# --------------------------------------------------------------------------
# LOD distance: pred_distance_lod() in appv6+.R
# --------------------------------------------------------------------------

def lod_distance(table: pd.DataFrame, maxf: pd.DataFrame, LOD: float) -> pd.DataFrame:
    """x_LOD = (ln LOD - intercept - lnmax) / slope, with avg using mean(lnmax),
    upper using (conf.low intercept, conf.high slope, max lnmax), lower using
    (conf.high intercept, conf.low slope, min lnmax). Negative -> NaN."""
    lod = np.log(LOD)
    rows = []
    for f in table["Filter"].unique():
        ti = table[(table["Filter"] == f) & (table["term"] == "(Intercept)")].iloc[0]
        ts = table[(table["Filter"] == f) & (table["term"] == "Distance")].iloc[0]
        d = (lod - ti["estimate"] - maxf["lnmax"].mean()) / ts["estimate"]
        du = (lod - ti["conf.low"] - maxf["lnmax"].max()) / ts["conf.high"]
        dl = (lod - ti["conf.high"] - maxf["lnmax"].min()) / ts["conf.low"]
        fix = lambda v: np.nan if (not np.isfinite(v) or v < 0) else v
        rows.append(dict(Filter=f, avg=fix(d), upper=fix(du), lower=fix(dl)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Reactive SMIM: k_func(), C_func(), optim_param() in the original
# --------------------------------------------------------------------------

LOGT1, LOGT2 = 0.0, 5.0      # hard-coded in appv6+.R: t1 = 1 s, t2 = 1e5 s


def gammainc_upper(a: float, x) -> np.ndarray:
    """expint::gammainc(a, x) = upper incomplete gamma Gamma(a, x), valid for a < 0."""
    import mpmath
    x = np.atleast_1d(np.asarray(x, float))
    return np.array([float(mpmath.gammainc(a, xi, mpmath.inf)) for xi in x])


def k_func(rs, rh, U, D, Lambda, beta, logT1=LOGT1, logT2=LOGT2):
    """Eq. 4.15. Returns k > 0 with C = C0 exp(-k x)."""
    t1, t2 = 10 ** logT1, 10 ** logT2
    num = np.exp(rh * t1) * gammainc_upper(-beta, t1 * (1 / t2 + rh)) * (1 + t2 * rh) ** beta
    den = gammainc_upper(-beta, t1 / t2)
    Theta = rs + Lambda * (1 - num / den)
    k = (-U + np.sqrt(U ** 2 + 4 * D * Theta)) / (2 * D)
    k = np.asarray(k)
    return float(k.reshape(-1)[0]) if k.size == 1 else k


def fit_rh(k_data, rs, U, D, Lambda, beta, logT1=LOGT1, logT2=LOGT2, bounds=(1e-9, 1e2)):
    """optim_param(): minimise |k_func(rs, rh, ...) - k_data| over rh >= 0.
    The original uses optim(method='SANN') from rh = 1e-3; here a bounded
    deterministic search on log10(rh) is used (same objective, same answer when
    SANN converges). Returns (rh, residual)."""
    from scipy.optimize import minimize_scalar
    obj = lambda lr: abs(float(k_func(rs, 10 ** lr, U, D, Lambda, beta, logT1, logT2)) - k_data)
    r = minimize_scalar(obj, bounds=(np.log10(bounds[0]), np.log10(bounds[1])), method="bounded",
                        options={"xatol": 1e-6})
    return 10 ** r.x, r.fun
