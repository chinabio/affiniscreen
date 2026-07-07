#!/usr/bin/env python
"""validate_restraint_leg_inputs.py (v2.5.42)

Pre-submission validator for the Option A single-copy restraint leg. Scans
every mdin (lambda_*/*.in) in a leg directory and fails if any TI keyword is
present -- the restraint leg must be plain MD (icfe=0, lambda-scaled Boresch
&rst), so icfe=1 / timask / scmask / clambda / *mbar* / crgmask must NOT appear.

Catches the class of bugs seen in v2.5.36-2.5.41 where a single stage template
(heat) or the prod recovery wrapper still assumed TI -- before a GPU run is
wasted.

Usage:
    python tools/validate_restraint_leg_inputs.py <leg_dir> [<leg_dir> ...]

Exit 0 = all clean, 1 = TI keyword found, 2 = usage / no mdin found.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import sys
from pathlib import Path

FORBIDDEN = ("icfe=1", "timask1", "timask2", "scmask1", "scmask2",
             "clambda", "ifmbar", "mbar_lambda", "mbar_states", "crgmask")


def check_leg(leg_dir):
    mdins = sorted(Path(leg_dir).glob("lambda_*/*.in"))
    if not mdins:
        print(f"  WARNING: no lambda_*/*.in files under {leg_dir}")
        return None
    offenders = []
    for mdin in mdins:
        try:
            text = mdin.read_text()
        except OSError:
            continue
        hits = [kw for kw in FORBIDDEN if kw in text]
        if hits:
            offenders.append((mdin, hits))
    if offenders:
        print(f"  FAIL: {leg_dir}")
        for mdin, hits in offenders:
            print(f"    {mdin.relative_to(leg_dir)}: {', '.join(hits)}")
        return False
    print(f"  OK ({len(mdins)} mdin clean): {leg_dir}")
    return True


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    results = [check_leg(d) for d in argv[1:]]
    if any(r is False for r in results):
        print("\nRESULT: TI keywords found in a restraint leg -- do NOT submit.")
        return 1
    if all(r is None for r in results):
        return 2
    print("\nRESULT: all restraint-leg inputs clean (Option A).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
