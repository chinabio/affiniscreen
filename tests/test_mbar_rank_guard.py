"""v2.5.15 regression: MBAR rank guard."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import numpy as np
import pandas as pd
from amber_md.fep import FEPAnalyzer, _UNK_COND_MAX


def _mk(cols, rows=40, deficient=False, seed=0):
    rng = np.random.default_rng(seed)
    st = list(cols); blocks = []
    for lam in st:
        b = np.array([[abs(lam - mu) * 6.0 for mu in st] for _ in range(rows)])
        b = b + rng.normal(0, 0.05, b.shape)
        if deficient:
            b[:, 2:] = b[:, :1]
        idx = pd.MultiIndex.from_tuples([(float(t), lam) for t in range(rows)],
                                        names=["time", "lambdas"])
        blocks.append(pd.DataFrame(b, columns=st, index=idx))
    return pd.concat(blocks)


def test_full_rank_well_posed():
    cols = [0.0, 0.25, 0.5, 0.75, 1.0]
    r, c, n, well = FEPAnalyzer._u_nk_rank_cond(_mk(cols))
    assert n == len(cols) and r == len(cols) and c <= _UNK_COND_MAX and well


def test_rank_deficient_flagged():
    cols = [0.0, 0.25, 0.5, 0.75, 1.0]
    r, c, n, well = FEPAnalyzer._u_nk_rank_cond(_mk(cols, deficient=True))
    assert r < n and well is False


def test_empty_never_raises():
    r, c, n, well = FEPAnalyzer._u_nk_rank_cond(pd.DataFrame())
    assert well is True
