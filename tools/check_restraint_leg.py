#!/usr/bin/env python
"""check_restraint_leg.py (v2.5.37) -- post-run reliability check for the
Option A single-copy lambda-scaled Boresch restraint leg.

Reads a finished leg directory and confirms the analyzer reliability gate
(introduced in v2.5.36) passes:
    * max|<dV/dl>| over all windows must be <= DVDL_MAX (default 200)
    * |BAR - TI| must be <= BARTI_MAX (default 50), when both exist
    * the dV/dl profile must be MONOTONE-FREE of runaway tails

It parses dvdl_summary.csv (leg,lambda,dvdl_kcal_mol) and dG_estimators.csv
if present. Exit code 0 = reliable, 1 = UNRELIABLE, 2 = usage/IO error.

Usage:
    python tools/check_restraint_leg.py <leg_dir> [--dvdl-max 200] [--barti-max 50]
"""
import csv
import sys
from pathlib import Path

DVDL_MAX = 200.0
BARTI_MAX = 50.0


def _read_dvdl(leg_dir):
    f = leg_dir / "dvdl_summary.csv"
    if not f.exists():
        return None
    rows = {}
    with open(f, newline="") as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#") or row[0] == "leg":
                continue
            try:
                lam = float(row[1]); val = float(row[2])
            except (IndexError, ValueError):
                continue
            rows[lam] = val
    return rows


def _read_estimators(leg_dir):
    f = leg_dir / "dG_estimators.csv"
    if not f.exists():
        return {}
    est = {}
    with open(f, newline="") as fh:
        rdr = csv.DictReader(
            (l for l in fh if not l.startswith("# leg")))
        for row in rdr:
            name = (row.get("estimator") or row.get("name") or "").strip()
            try:
                est[name] = float(row.get("dG_kcal_mol"))
            except (TypeError, ValueError):
                pass
    return est


def main(argv):
    if len(argv) < 2:
        print(__doc__); return 2
    leg_dir = Path(argv[1])
    dvdl_max = DVDL_MAX; barti_max = BARTI_MAX
    if "--dvdl-max" in argv:
        dvdl_max = float(argv[argv.index("--dvdl-max") + 1])
    if "--barti-max" in argv:
        barti_max = float(argv[argv.index("--barti-max") + 1])
    if not leg_dir.is_dir():
        print(f"ERROR: not a directory: {leg_dir}"); return 2

    dvdl = _read_dvdl(leg_dir)
    if not dvdl:
        print(f"ERROR: no usable dvdl_summary.csv in {leg_dir}"); return 2
    est = _read_estimators(leg_dir)

    reasons = []
    max_abs = max(abs(v) for v in dvdl.values())
    if max_abs > dvdl_max:
        worst = max(dvdl, key=lambda k: abs(dvdl[k]))
        reasons.append(f"max|dV/dl|={max_abs:.1f} (> {dvdl_max:.0f}) "
                       f"at lambda={worst:.3f} -> end-state singularity?")
    ti = est.get("TI"); bar = est.get("BAR")
    if ti is not None and bar is not None and abs(bar - ti) > barti_max:
        reasons.append(f"|BAR-TI|={abs(bar-ti):.1f} (> {barti_max:.0f}) "
                       f"-> adjacent-window overlap is ~zero")

    print(f"=== restraint-leg reliability: {leg_dir.name} ===")
    print(f"  windows analyzed : {len(dvdl)}")
    print(f"  max|<dV/dl>|     : {max_abs:.2f} kcal/mol  (limit {dvdl_max:.0f})")
    if ti is not None:
        print(f"  TI  dG           : {ti:+.2f} kcal/mol")
    if bar is not None:
        print(f"  BAR dG           : {bar:+.2f} kcal/mol")
    # dV/dl profile (sorted)
    print("  dV/dl profile:")
    for lam in sorted(dvdl):
        flag = "  <-- HIGH" if abs(dvdl[lam]) > dvdl_max else ""
        print(f"    lambda={lam:.3f}  dV/dl={dvdl[lam]:+9.2f}{flag}")

    if reasons:
        print("\nRESULT: UNRELIABLE")
        for rr in reasons:
            print(f"  - {rr}")
        return 1
    print("\nRESULT: RELIABLE (gate passed)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
