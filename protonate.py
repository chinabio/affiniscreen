#!/usr/bin/env python
"""Preview / apply protonation state assignments."""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import sys, os, argparse
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amber_md.protonation import detect_protonation, apply_protonation

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pdb", type=Path)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--override", action="append", default=[])
    a = p.parse_args()
    print(f"Scanning {a.pdb} ...")
    suggestions = detect_protonation(a.pdb)
    if not suggestions:
        print("  No protonation changes detected.")
    else:
        print(f"  {len(suggestions)} residue(s) need renaming:")
        print(f"  {'CHAIN':<7}{'RES':<6}{'OLD':<5}{'NEW':<5}{'REASON'}")
        for c, r, old, new, why in suggestions:
            print(f"  {c or '_':<7}{r:<6}{old:<5}{new:<5}{why}")
    if a.output:
        overrides = {}
        for spec in a.override:
            chain, resnum, name = spec.split(":")
            overrides[(chain, int(resnum))] = name.upper()
        apply_protonation(a.pdb, a.output, manual_overrides=overrides or None)
        print(f"\nWrote: {a.output}")

if __name__ == "__main__":
    main()
