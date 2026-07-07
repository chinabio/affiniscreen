"""Regression tests for the final60 FEP estimator + eq-stability fixes.

These tests use REAL pmemd `eq.out` fixtures captured from the failing
`complex_vdw` leg of run abfe_20260609_122733:

    fixtures/eq_complex_vdw_lambda_0.450.out  -- healthy eq (peak ~300.7 K)
    fixtures/eq_complex_vdw_lambda_0.500.out  -- detonated eq (~14,963 K;
                                                 MBAR block is all '****')

They lock in two behaviours:

  BUG 2  -- the original `_sanitize_u_nk` dropped whole frames, collapsing the
            reduced-potential matrix and raising
            "Shape of passed values is (n,n), indices imply (n+1,n+1)".
            The patched clamp + `_reconcile_u_nk_states` must keep the matrix
            square so MBAR/BAR fit.

  BUG 1  -- a soft-core blow-up integrates to a non-physical temperature while
            still exiting 0. The eq stability gate must flag the λ=0.50 file
            and pass the λ=0.45 file, using the same peak-TEMP parse the
            generated shell gate uses.

Run:  pytest -q test_fep_estimator_recovery.py
"""
import os
import re
import numpy as np
import pandas as pd
import pytest

from alchemlyb.parsing.amber import extract_u_nk
from alchemlyb.estimators import MBAR, BAR
from alchemlyb.preprocessing import statistical_inefficiency

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
EQ_GOOD = os.path.join(FIX, "eq_complex_vdw_lambda_0.450.out")
EQ_BLOWN = os.path.join(FIX, "eq_complex_vdw_lambda_0.500.out")
T = 298.0
_UNK_SANE_MAX_KT = 1.0e4
EQ_TEMP_MAX_K = 1000.0


# --------------------------------------------------------------------------- #
# Reference implementations mirroring amber_md/fep.py (final60).               #
# Kept inline so the test is runnable standalone without importing the whole   #
# workflow package; they are byte-for-byte equivalent to the patched methods.  #
# --------------------------------------------------------------------------- #
def sanitize_u_nk_OLD(u_nk):
    """The buggy pre-final60 behaviour: DROP offending rows."""
    vals = u_nk.to_numpy(dtype=float, copy=True)
    if vals.size == 0:
        return u_nk, 0
    baseline = np.nanmedian(vals, axis=1, keepdims=True)
    baseline = np.where(np.isfinite(baseline), baseline, 0.0)
    bad = (~np.isfinite(vals)).any(axis=1) | (
        np.abs(vals - baseline) > _UNK_SANE_MAX_KT).any(axis=1)
    n = int(bad.sum())
    return (u_nk if n == 0 else u_nk.loc[~bad]), n


def sanitize_u_nk_NEW(u_nk):
    """final60: CLAMP in place, never drop rows."""
    vals = u_nk.to_numpy(dtype=float, copy=True)
    if vals.size == 0:
        return u_nk, 0
    baseline = np.nanmedian(vals, axis=1, keepdims=True)
    gmed = np.nanmedian(vals[np.isfinite(vals)]) if np.isfinite(vals).any() else 0.0
    baseline = np.where(np.isfinite(baseline), baseline, gmed)
    lo, hi = baseline - _UNK_SANE_MAX_KT, baseline + _UNK_SANE_MAX_KT
    nonfinite = ~np.isfinite(vals)
    over = np.isfinite(vals) & ((vals > hi) | (vals < lo))
    n = int(nonfinite.sum() + over.sum())
    if n:
        vals = np.where(nonfinite, np.broadcast_to(baseline, vals.shape), vals)
        vals = np.clip(vals, lo, hi)
        u_nk = u_nk.copy()
        u_nk.iloc[:, :] = vals
    return u_nk, n


def reconcile_u_nk_states(u_nk):
    """final60: trim rows+cols to the intersection so the matrix is square."""
    sampled = set(dict.fromkeys(u_nk.index.get_level_values("lambdas")))

    def _key(c):
        return c[0] if isinstance(c, tuple) and len(c) == 1 else c

    cols = list(u_nk.columns)
    keep = [c for c in cols if _key(c) in sampled]
    dropped = [c for c in cols if _key(c) not in sampled]
    if not dropped:
        return u_nk, []
    keep_keys = set(_key(c) for c in keep)
    u_nk = u_nk[keep]
    mask = u_nk.index.get_level_values("lambdas").isin(keep_keys)
    return u_nk.loc[mask], dropped


def decorrelate_per_lambda(df):
    out = []
    for _lam, g in df.groupby(level="lambdas"):
        g = g[~g.index.duplicated(keep="first")].sort_index()
        try:
            g = statistical_inefficiency(g, g.iloc[:, 0])
        except Exception:
            pass
        out.append(g)
    return pd.concat(out) if out else df


def peak_temp_from_mdout(path):
    """Mirror of the generated shell gate's TEMP(K) parse.

    Matches ONLY the main-system energy lines (those beginning with
    ' NSTEP ='). pmemd softcore runs also print a per-step
    'Softcore part of the system:  N atoms,  TEMP(K) = ...' line whose
    small-subsystem temperature is noisy (routinely 200-410 K and can
    transiently exceed 1000 K on a perfectly healthy window); including
    those would cause false-positive eq failures. The shell gate uses
    the equivalent `grep -aE '^ *NSTEP ='` pre-filter.
    """
    temps = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.lstrip().startswith("NSTEP ="):
                continue
            for m in re.findall(r"TEMP\(K\) *= *([-0-9.]+)", line):
                try:
                    temps.append(float(m))
                except ValueError:
                    pass
    return max(temps) if temps else None


# --------------------------------------------------------------------------- #
# Fixtures                                                                     #
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def states():
    u = extract_u_nk(EQ_BLOWN, T=T)
    return list(u.columns)


def _healthy_window(lam, states, M=80, seed=0):
    rng = np.random.default_rng(seed)
    rows = [[abs(lam - mu) * 6.0 + rng.normal(0, 0.3) for mu in states]
            for _ in range(M)]
    idx = pd.MultiIndex.from_tuples(
        [(float(t), lam) for t in range(M)], names=["time", "lambdas"])
    return pd.DataFrame(rows, columns=states, index=idx)


def _assemble_leg(states, corrupt_fraction=0.95, seed=1):
    """A leg: 3 healthy windows + 1 mostly-corrupt window (the real λ=0.50)."""
    u050 = extract_u_nk(EQ_BLOWN, T=T)
    arr = u050.to_numpy(float).copy()
    rng = np.random.default_rng(seed)
    kill = rng.random(arr.shape[0]) < corrupt_fraction
    arr[kill, :] = np.inf
    bad = pd.DataFrame(arr, index=u050.index, columns=u050.columns)
    parts = [_healthy_window(l, states, seed=i)
             for i, l in enumerate((0.0, 0.05, 0.45))]
    parts.append(bad)
    return pd.concat(parts)


# --------------------------------------------------------------------------- #
# BUG 2: estimator-matrix recovery                                            #
# --------------------------------------------------------------------------- #
def test_old_sanitizer_reproduces_shape_error(states):
    """The pre-final60 drop-rows path must still reproduce the (n,n)/(n+1,n+1)
    failure -- this guards against accidentally 'fixing' the bug in the OLD
    reference and silently weakening the regression."""
    u = _assemble_leg(states)
    u_old, _ = sanitize_u_nk_OLD(u)
    u_old = decorrelate_per_lambda(u_old)
    with pytest.raises(ValueError, match=r"indices imply"):
        MBAR().fit(u_old)


def test_patched_pipeline_recovers_mbar_and_bar(states):
    """final60: clamp + reconcile yields a square matrix and BOTH estimators
    fit on real corrupt data."""
    u = _assemble_leg(states)
    u_new, n_clamped = sanitize_u_nk_NEW(u)
    assert n_clamped > 0                         # corruption was present
    assert u_new.shape[0] == u.shape[0]          # NO rows dropped
    u_new = decorrelate_per_lambda(u_new)
    u_new, dropped = reconcile_u_nk_states(u_new)

    # square: every column has a matching sampled lambda group
    sampled = set(u_new.index.get_level_values("lambdas"))
    assert set(u_new.columns) == sampled

    kT = 0.001987204 * T
    for Est in (BAR, MBAR):
        est = Est().fit(u_new)                   # must not raise
        dG = float(est.delta_f_.iloc[0, -1] * kT)
        assert np.isfinite(dG)


def test_clamp_preserves_all_finite_when_clean(states):
    """A clean window must pass through untouched (n_clamped == 0)."""
    clean = _healthy_window(0.30, states, M=50)
    out, n = sanitize_u_nk_NEW(clean)
    assert n == 0
    pd.testing.assert_frame_equal(out, clean)


# --------------------------------------------------------------------------- #
# BUG 1: eq stability gate                                                     #
# --------------------------------------------------------------------------- #
def test_blown_window_flagged_by_gate():
    peak = peak_temp_from_mdout(EQ_BLOWN)
    assert peak is not None
    assert peak > EQ_TEMP_MAX_K          # ~14,963 K -> gate fires
    assert peak > 10000.0                # sanity: this really is the detonation


def test_healthy_window_passes_gate():
    peak = peak_temp_from_mdout(EQ_GOOD)
    assert peak is not None
    assert peak < EQ_TEMP_MAX_K          # ~300.7 K -> gate does NOT fire
    assert 250.0 < peak < 320.0          # sanity: equilibrated at setpoint


def test_gate_threshold_separates_the_two():
    """The default threshold must sit strictly between the two real windows."""
    good = peak_temp_from_mdout(EQ_GOOD)
    blown = peak_temp_from_mdout(EQ_BLOWN)
    assert good < EQ_TEMP_MAX_K < blown