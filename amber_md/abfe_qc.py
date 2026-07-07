#!/usr/bin/env python3
"""abfe_qc.py -- QC report for an Amber ABFE work-dir.

Reads, for a finished (or partial) ABFE run:
  <wd>/fep/ABFE_RESULT.json          (cycle-closer output)
  <wd>/fep/ABFE_RESULT.txt           (human summary, if present)
  <wd>/fep/<leg>/summary.json        (per-leg analyzer: dG, overlap, windows)
where <leg> in: complex_decharge complex_vdw solvent_decharge solvent_vdw

Reports: cycle closure, per-leg dG breakdown, MBAR overlap adequacy,
window completeness, estimator spread / convergence, and a GO/REVIEW verdict.

Usage:
    python abfe_qc.py <work-dir>          # e.g. ~/abfe_smoketest_2026...
    python abfe_qc.py <work-dir> --json   # machine-readable
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import argparse, json, sys
from pathlib import Path

LEGS = ["complex_decharge", "complex_vdw", "solvent_decharge", "solvent_vdw"]
OVERLAP_WARN = 0.03   # min adjacent-window MBAR overlap considered adequate

def _load(p):
    try: return json.loads(Path(p).read_text())
    except Exception: return None

def _fmt(x, u="kcal/mol"):
    return "n/a" if x is None else f"{float(x):+.2f} {u}"

def main(argv=None):
    ap = argparse.ArgumentParser(prog="abfe_qc")
    ap.add_argument("work_dir", type=Path)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)

    fep = Path(a.work_dir).expanduser() / "fep"
    report = {"work_dir": str(a.work_dir), "checks": [], "verdict": "GO"}
    def add(level, msg):
        report["checks"].append({"level": level, "msg": msg})
        if level == "BLOCK": report["verdict"] = "REVIEW"

    if not fep.is_dir():
        print(f"ERROR: {fep} not found. Has the run started?"); return 2

    # ---- 1. cycle-closer result ----
    res = _load(fep / "ABFE_RESULT.json")
    if res is None:
        add("BLOCK", "ABFE_RESULT.json missing -- cycle closer has not written a "
                     "result yet (job still running, or a leg failed).")
    else:
        report["dG_bind_kcal_mol"] = res.get("dG_bind_kcal_mol")
        report["temperature_K"] = res.get("temperature_K")
        if res.get("trusted") is True:
            add("OK", f"Cycle closed & TRUSTED. dG_bind = "
                      f"{_fmt(res.get('dG_bind_kcal_mol'))}")
        else:
            add("BLOCK", f"Cycle closed but UNTRUSTED. dG_bind = "
                         f"{_fmt(res.get('dG_bind_kcal_mol'))}; "
                         f"incomplete={res.get('incomplete_legs')}, "
                         f"nonconverged={res.get('nonconverged_legs')}")
        add("INFO", f"complex+restr = {_fmt(res.get('dG_complex_plus_restr_kcal_mol'))} "
                    f"| solvent = {_fmt(res.get('dG_solvent_kcal_mol'))}")

    # ---- 2. per-leg analyzer summaries ----
    report["legs"] = {}
    for leg in LEGS:
        s = _load(fep / leg / "summary.json")
        if s is None:
            add("BLOCK", f"[{leg}] summary.json missing (leg not analyzed / failed).")
            continue
        nw, nr = s.get("n_windows"), s.get("n_requested")
        dG = s.get("dG_kcal_mol")
        est = s.get("headline_estimator") or s.get("estimator_used")
        ov = s.get("overlap")
        report["legs"][leg] = {"dG_kcal_mol": dG, "n_windows": nw,
                               "n_requested": nr, "estimator": est, "overlap": ov}
        # window completeness
        if nr and nw is not None and nw < nr:
            add("BLOCK", f"[{leg}] only {nw}/{nr} windows analyzed "
                         f"(missing: {s.get('missing_windows')}).")
        else:
            add("OK", f"[{leg}] dG={_fmt(dG)} via {est}; windows {nw}/{nr}.")
        # MBAR overlap adequacy
        worst = _min_overlap(ov)
        if est and "mbar" in str(est).lower():
            if worst is None:
                add("WARN", f"[{leg}] MBAR used but no overlap matrix found.")
            elif worst < OVERLAP_WARN:
                add("WARN", f"[{leg}] weak MBAR overlap (min adjacent "
                            f"{worst:.3f} < {OVERLAP_WARN}); add lambda windows here.")
            else:
                add("OK", f"[{leg}] MBAR overlap OK (min adjacent {worst:.3f}).")
        elif est:
            add("INFO", f"[{leg}] estimator '{est}' (no MBAR overlap diagnostic).")
        # convergence / estimator spread (final40: distinguish a genuine
        # TI-vs-MBAR disagreement from an outright MBAR/BAR FAILURE).
        if s.get("only_ti_survived"):
            ff = s.get("failed_estimators") or {}
            why = "; ".join(f"{k}: {v}" for k, v in ff.items()) or "SVD/parse failure"
            add("WARN", f"[{leg}] MBAR/BAR unavailable -- TI-only result "
                        f"({why}). Typical of truncated or too-few-window runs; "
                        f"not a TI/MBAR disagreement.")
        else:
            sp = s.get("estimator_spread_kcal")
            if sp is not None and abs(float(sp)) > 1.0:
                add("WARN", f"[{leg}] estimator spread {float(sp):.2f} kcal/mol "
                            f"(>1.0) -- TI/MBAR disagree; lengthen production / "
                            f"add windows.")

    _print(report) if not a.json else print(json.dumps(report, indent=2, default=str))
    return 0 if report["verdict"] == "GO" else 1

def _min_overlap(ov):
    """Smallest adjacent-window overlap from a scalar or NxN matrix."""
    if ov is None: return None
    if isinstance(ov, (int, float)): return float(ov)
    try:
        m = ov
        if isinstance(m, dict): m = m.get("matrix", None)
        if not m: return None
        adj = [m[i][i+1] for i in range(len(m)-1)] + [m[i+1][i] for i in range(len(m)-1)]
        return min(float(x) for x in adj) if adj else None
    except Exception:
        return None

def _print(r):
    print("=" * 66)
    print(f"ABFE QC report  --  {r['work_dir']}")
    print("=" * 66)
    if r.get("dG_bind_kcal_mol") is not None:
        print(f"  dG_bind = {r['dG_bind_kcal_mol']:+.2f} kcal/mol"
              f"   (T={r.get('temperature_K')} K)\n")
    order = {"BLOCK": 0, "WARN": 1, "OK": 2, "INFO": 3}
    for c in sorted(r["checks"], key=lambda c: order.get(c["level"], 9)):
        tag = {"OK": "[ OK ]", "WARN": "[WARN]", "BLOCK": "[BLOCK]",
               "INFO": "[info]"}[c["level"]]
        print(f"  {tag:8s} {c['msg']}")
    print("\n" + "=" * 66)
    print(f"VERDICT: {r['verdict']}"
          + ("  -- result looks usable." if r["verdict"] == "GO"
             else "  -- inspect the [BLOCK]/[WARN] items above."))
    print("=" * 66)

if __name__ == "__main__":
    sys.exit(main())
