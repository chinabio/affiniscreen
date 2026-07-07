"""v2.5.15 regression: completeness gate + recovery logging."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import json
from pathlib import Path
from amber_md.recovery import (
    leg_completeness, ensure_leg_complete, write_recovery_report)

_DONE = "Master Total wall time   12.3 seconds\nFinal Performance Info\n"


def _win(leg, lam, complete):
    wd = Path(leg) / ("lambda_%.3f" % lam)
    wd.mkdir(parents=True, exist_ok=True)
    wd.joinpath("prod.out").write_text(_DONE if complete else "vlimit exceeded\n")
    return wd


def test_completeness_counts(tmp_path):
    leg = tmp_path / "complex_vdw"
    _win(leg, 0.0, True); _win(leg, 0.5, False); _win(leg, 1.0, True)
    c, n, missing = leg_completeness(leg, [0.0, 0.5, 1.0])
    assert (c, n) == (2, 3) and missing == ["0.500"]


def test_gate_recovers_and_logs(tmp_path):
    leg = tmp_path / "complex_vdw"
    _win(leg, 0.0, True); bad = _win(leg, 0.5, False); _win(leg, 1.0, True)
    # stub rerun: "completes" the bad window by writing a done marker.
    calls = []
    def rerun(window_dir, attempt):
        calls.append((Path(window_dir).name, attempt))
        (Path(window_dir) / "prod.out").write_text(_DONE)
        return 0
    rec = ensure_leg_complete(leg, [0.0, 0.5, 1.0], rerun, max_attempts=3)
    assert rec["recovery_attempted"] is True
    assert rec["pre"]["missing"] == ["0.500"]
    assert rec["post"]["missing"] == []
    assert rec["complete"] is True
    # structured log written + parseable
    log = json.loads((leg / "recovery_log.json").read_text())
    assert log["leg"] == "complex_vdw"
    assert "0.500" in log["windows"]
    assert log["windows"]["0.500"]["completed"] is True
    assert calls, "rerun was never invoked"


def test_gate_skips_when_already_complete(tmp_path):
    leg = tmp_path / "solvent"
    _win(leg, 0.0, True); _win(leg, 1.0, True)
    def rerun(window_dir, attempt):
        raise AssertionError("rerun must NOT be called when all complete")
    rec = ensure_leg_complete(leg, [0.0, 1.0], rerun, max_attempts=3)
    assert rec["recovery_attempted"] is False and rec["complete"] is True


def test_report_rollup(tmp_path):
    leg = tmp_path / "complex_vdw"
    _win(leg, 0.0, True)
    rec = ensure_leg_complete(leg, [0.0], lambda *a: 0, max_attempts=1)
    ok = write_recovery_report(tmp_path, [rec])
    txt = (tmp_path / "RECOVERY_REPORT.txt").read_text()
    assert "AMBER FEP recovery report" in txt and "LEG: complex_vdw" in txt
    assert ok is True
