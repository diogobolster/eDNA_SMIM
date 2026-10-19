"""Checks of edna_py/original.py, the Python port of the R application's model code.

`test_r_*` compare against numbers computed in R (tests/r_reference.tsv: lm() with
confint(), expint::gammainc(), Eq. 4.15 and its inversion on the application's template data).
They are what shows the port is faithful rather than merely self-consistent. The remaining
tests are properties the implementation has to satisfy on its own terms, including the
conventions of the original code that the thesis text leaves open.

Run from the repository root:  python -m pytest tests
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "edna_py"))
import original as O  # noqa: E402

REL = 1e-6


@pytest.fixture(scope="module")
def ref() -> dict:
    out = {}
    for line in (ROOT / "tests" / "r_reference.tsv").read_text().splitlines():
        if "\t" in line:
            k, v = line.split("\t")
            out[k] = float(v)
    return out


@pytest.fixture(scope="module")
def template() -> pd.DataFrame:
    return pd.read_csv(ROOT / "app_R" / "template_sampling_perSize.csv", encoding="utf-8-sig")


@pytest.fixture(scope="module")
def source() -> pd.DataFrame:
    return pd.read_csv(ROOT / "app_R" / "template_source.csv", encoding="utf-8-sig")


@pytest.fixture(scope="module")
def sampling(template) -> pd.DataFrame:
    return template[(template["Treatment"] == "A_NP_STHD") & (template["Target"] == "STHD")]


# --- regression against R's lm() + confint() --------------------------------

def test_r_regression_per_component(sampling, ref):
    for filt in (0.2, 1.0, 10.0):
        s = sampling[(sampling["Filter"] == filt) & (sampling["Concentration"] > 0)]
        fit = O._lm(s["Distance"], np.log(s["Concentration"]))
        pre = f"fit.{filt:g}."
        slope, icpt = fit.iloc[1], fit.iloc[0]
        assert -slope["estimate"] == pytest.approx(ref[pre + "lambda"], rel=REL)
        assert -slope["conf.high"] == pytest.approx(ref[pre + "lambda_low"], rel=REL)
        assert -slope["conf.low"] == pytest.approx(ref[pre + "lambda_up"], rel=REL)
        assert np.exp(icpt["estimate"]) == pytest.approx(ref[pre + "C0"], rel=REL)
        assert np.exp(icpt["conf.low"]) == pytest.approx(ref[pre + "C0_low"], rel=REL)
        assert np.exp(icpt["conf.high"]) == pytest.approx(ref[pre + "C0_up"], rel=REL)
        assert slope["r.squared"] == pytest.approx(ref[pre + "r2"], rel=REL)
        assert slope["n"] == ref[pre + "n"]


def test_r_regression_total_and_lod(sampling, ref):
    tot = sampling.groupby(["Replicate", "Distance"], as_index=False)["Concentration"].sum()
    tot = tot[tot["Concentration"] > 0]
    fit = O._lm(tot["Distance"], np.log(tot["Concentration"]))
    slope, icpt = fit.iloc[1], fit.iloc[0]
    assert -slope["estimate"] == pytest.approx(ref["total.lambda"], rel=REL)
    assert np.exp(icpt["estimate"]) == pytest.approx(ref["total.C0"], rel=REL)
    assert slope["r.squared"] == pytest.approx(ref["total.r2"], rel=REL)
    # LOD distance x = ln(C0 / LOD) / lambda with LOD = 0.05, best estimate and the two
    # confidence-interval corners (low C0 with high lambda, high C0 with low lambda)
    LOD = 0.05
    best = (icpt["estimate"] - np.log(LOD)) / -slope["estimate"]
    low = (icpt["conf.low"] - np.log(LOD)) / -slope["conf.low"]
    up = (icpt["conf.high"] - np.log(LOD)) / -slope["conf.high"]
    assert best == pytest.approx(ref["lod.best"], rel=REL)
    assert low == pytest.approx(ref["lod.low"], rel=REL)
    assert up == pytest.approx(ref["lod.up"], rel=REL)


# --- SMIM against R (expint::gammainc, Eq. 4.15, inversion) --------------------

V, D, LAM, BETA = 0.0982, 0.0345, 0.085, 0.898
LOGT1, LOGT2 = 0.0, 3.0          # t1 = 1 s, t2 = 1000 s in the R reference run


def test_r_upper_incomplete_gamma_negative_order(ref):
    for x in (0.001, 0.01, 0.1, 0.5, 1, 5):
        assert O.gammainc_upper(-0.898, x)[0] == pytest.approx(ref[f"gamma_upper.-0.898.{x:g}"], rel=REL)
    for x in (0.2, 2.0):
        assert O.gammainc_upper(-1.5, x)[0] == pytest.approx(ref[f"gamma_upper.-1.5.{x:g}"], rel=REL)


def test_r_k_of_eq_4_15(ref):
    # R reference stores k with the sign convention C = C0 exp(k x), k < 0; the port returns
    # k > 0 with C = C0 exp(-k x).
    for rh in (1e-5, 1e-3, 1e-1):
        k = O.k_func(3e-5, rh, V, D, LAM, BETA, LOGT1, LOGT2)
        assert k == pytest.approx(-ref[f"k.{rh:g}"], rel=REL)
    # Theta = 0 when rs = rh = 0: k = 0 exactly
    assert O.k_func(0.0, 0.0, V, D, LAM, BETA, LOGT1, LOGT2) == pytest.approx(0.0, abs=1e-15)


def test_r_rh_inversion_round_trip(ref):
    k = O.k_func(3e-5, 1e-3, V, D, LAM, BETA, LOGT1, LOGT2)
    rh, resid = O.fit_rh(k, 3e-5, V, D, LAM, BETA, LOGT1, LOGT2)
    assert rh == pytest.approx(ref["fit_rh"], rel=1e-4)
    assert rh == pytest.approx(1e-3, rel=1e-4)
    assert resid < 1e-7


def test_k_increases_with_rh():
    ks = [O.k_func(2e-5, rh, V, D, LAM, BETA) for rh in (1e-6, 1e-4, 1e-2, 1.0)]
    assert all(np.diff(ks) > 0)


# --- conventions of the original application ------------------------------------

def test_normalise_replicate_max_and_total_filter(template):
    d = O.normalise(template)
    assert set(["max_c", "norm_c", "ln_norm_c"]) <= set(d.columns)
    assert (d["norm_c"] <= 1 + 1e-12).all()
    # each replicate's maximum normalises to 1
    m = d.groupby(["Treatment", "Replicate"])["norm_c"].max()
    assert np.allclose(m, 1.0)
    # total data without a Filter column gets one
    t = O.normalise(template.drop(columns=["Filter"]).groupby(
        ["Treatment", "Replicate", "Target", "Distance"], as_index=False)["Concentration"].sum())
    assert (t["Filter"] == "total").all()


def test_source_composition_sums_to_one(source):
    f = O.source_composition(source[source["Target"] == "STHD"])
    assert f["mean"].sum() == pytest.approx(1.0)
    assert (f["std"] >= 0).all()


@pytest.fixture(scope="module")
def mc(template):
    d = O.normalise(template)
    tab = O.table_var(d, "A_NP_STHD", by_size=True)
    return d, tab, O.mc_new(tab, nsim=400, dt=0.5, tmax=100, seed=1)


def test_mc_slope_is_clipped_at_zero(mc):
    _, _, m = mc
    # a drawn slope that is positive (concentration rising downstream) is replaced by zero
    assert (m["slope_used"] <= 0).all()
    assert np.array_equal(m["slope_used"], np.minimum(m["slope"], 0))
    # fractions over the components (axis 1) sum to one at every grid point
    assert m["frac"].shape[1] == len(m["filters"])
    assert np.allclose(m["frac"].sum(axis=1), 1.0)


def test_mc_solutions_respect_tolerance_and_window(mc, source):
    d, tab, m = mc
    f0 = O.source_composition(source[source["Target"] == "STHD"])["mean"].to_numpy()
    sol = O.df_close(m, f0, tol=0.05, tmax=100, denorm=O.denorm_factor(d, "A_NP_STHD"))
    assert len(sol) > 0
    assert (sol["distance"] < 100).all() and (sol["distance"] > 0).all()     # grid ends discarded
    assert (sol["SSe"] < 0.05).all()
    assert np.allclose(sol["C0denorm"], sol["C0"] * O.denorm_factor(d, "A_NP_STHD"))


def test_mc_tighter_tolerance_accepts_fewer(mc, source):
    d, tab, m = mc
    f0 = O.source_composition(source[source["Target"] == "STHD"])["mean"].to_numpy()
    n = [len(O.df_close(m, f0, tol=t, tmax=100)) for t in (0.2, 0.05, 0.01)]
    assert n[0] >= n[1] >= n[2]


def test_mc_seed_is_reproducible(mc):
    _, tab, m1 = mc
    m2 = O.mc_new(tab, nsim=400, dt=0.5, tmax=100, seed=1)
    assert np.array_equal(m1["frac"], m2["frac"])


def test_hist_mode_is_bin_centre():
    x = np.concatenate([np.full(50, 3.0), np.linspace(0, 15, 50)])
    mode = O.hist_mode(x, bins=15)
    assert abs(mode - 3.0) <= 0.5


def test_two_components_make_tolerance_irrelevant(template, source):
    """With n = 2 the composition is one number and a trajectory is monotonic in it, so its
    closest approach to the source composition is the grid point nearest the crossing (or a
    window end, discarded). The tolerance therefore plays no part once it exceeds half the
    composition change over one grid step; a smaller tolerance only removes trajectories whose
    crossing falls between grid points, it never moves a solution."""
    t = template[(template["Treatment"] == "A_NP_STHD") & (template["Filter"] != 10.0)]
    d = O.normalise(t)
    tab = O.table_var(d, "A_NP_STHD", by_size=True)
    m = O.mc_new(tab, nsim=300, dt=0.5, tmax=100, seed=3)
    s = source[(source["Target"] == "STHD") & (source["Filter"] != 10.0)]
    f0 = O.source_composition(s)["mean"].to_numpy()
    half_step = 0.5 * np.abs(np.diff(m["frac"][:, 0, :], axis=-1)).max()
    assert half_step < 0.05
    a = O.df_close(m, f0, tol=0.05, tmax=100)
    b = O.df_close(m, f0, tol=1.0, tmax=100)
    assert len(a) == len(b) > 0
    assert np.allclose(a["distance"].to_numpy(), b["distance"].to_numpy())
    c = O.df_close(m, f0, tol=0.005, tmax=100)          # below the grid resolution
    assert len(c) < len(a)
    j = a.set_index("Nsim").join(c.set_index("Nsim"), rsuffix="_c", how="inner")
    assert len(j) == len(c)
    assert np.allclose(j["distance"], j["distance_c"])


# --- the thesis validation, recomputed ----------------------------------------

def test_mesocosm_validation_results_are_present():
    """run_stage0b_mesocosm.py writes the recomputation of thesis Table 4.2; run A must bracket
    the known 5.33 h in at least 8 of 12 datasets (the thesis reports 9 of 12)."""
    p = ROOT / "results" / "stage0b_mesocosm_table42_recomputed.csv"
    if not p.exists():
        pytest.skip("run run_all.sh first")
    r = pd.read_csv(p)
    a = r[r["run"] == "A_thesis_f_eps"]
    assert len(a) == 12
    cov = ((a["x_P10"] <= 5.33) & (a["x_P90"] >= 5.33)).sum()
    assert cov >= 8
