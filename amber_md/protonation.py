"""Auto-detect protein protonation states from atom presence."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
from .logger import get_logger
log = get_logger()


def detect_protonation(pdb_path):
    residue_atoms = defaultdict(set)
    residue_meta = {}
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("ATOM"): continue
        try: resnum = int(line[22:26])
        except ValueError: continue
        chain = line[21:22]; resname = line[17:20].strip(); atom = line[12:16].strip()
        key = (chain, resnum)
        residue_atoms[key].add(atom)
        residue_meta[key] = resname

    suggestions = []
    for key, atoms in residue_atoms.items():
        chain, resnum = key; resname = residue_meta[key]
        if resname == "HIS":
            has_hd1 = "HD1" in atoms; has_he2 = "HE2" in atoms
            if has_hd1 and has_he2:
                suggestions.append((chain, resnum, resname, "HIP",
                                    "doubly protonated (HD1+HE2 present, charge +1)"))
            elif has_hd1:
                suggestions.append((chain, resnum, resname, "HID",
                                    "delta-protonated (HD1 only, neutral)"))
            elif has_he2:
                suggestions.append((chain, resnum, resname, "HIE",
                                    "epsilon-protonated (HE2 only, neutral)"))
        elif resname == "GLU":
            if "HE2" in atoms or "HE1" in atoms:
                which = "HE2" if "HE2" in atoms else "HE1"
                suggestions.append((chain, resnum, resname, "GLH",
                                    f"protonated COOH ({which} present, neutral)"))
        elif resname == "ASP":
            if "HD2" in atoms or "HD1" in atoms:
                which = "HD2" if "HD2" in atoms else "HD1"
                suggestions.append((chain, resnum, resname, "ASH",
                                    f"protonated COOH ({which} present, neutral)"))
        elif resname == "LYS":
            if {"HZ1","HZ2"}.issubset(atoms) and "HZ3" not in atoms:
                suggestions.append((chain, resnum, resname, "LYN",
                                    "deprotonated NH2 (HZ3 absent, neutral)"))
    return suggestions


def apply_protonation(in_pdb, out_pdb, *, auto=True, manual_overrides=None):
    suggestions = detect_protonation(in_pdb) if auto else []
    suggestion_map = {(c, r): (old, new, why) for c, r, old, new, why in suggestions}
    if manual_overrides:
        for (c, r), new_name in manual_overrides.items():
            existing = suggestion_map.get((c, r))
            old_name = existing[0] if existing else "?"
            suggestion_map[(c, r)] = (old_name, new_name, "manual override")

    new_lines = []; applied = []; applied_keys = set()
    for line in Path(in_pdb).read_text().splitlines():
        if line.startswith("ATOM"):
            try:
                resnum = int(line[22:26]); chain = line[21:22]; old = line[17:20].strip()
                key = (chain, resnum)
                if key in suggestion_map:
                    real_old, new, why = suggestion_map[key]
                    if old == real_old or real_old == "?":
                        line = line[:17] + f"{new:>3s}" + line[20:]
                        if key not in applied_keys:
                            applied.append((chain, resnum, old, new, why))
                            applied_keys.add(key)
            except ValueError: pass
        new_lines.append(line)
    Path(out_pdb).write_text("\n".join(new_lines) + "\n")
    if applied:
        log.info("Protonation: renamed %d residue(s):", len(applied))
        for c, r, old, new, why in applied:
            log.info("  chain %s  %s %d -> %s  (%s)", c or "_", old, r, new, why)
    else:
        log.info("Protonation: no changes needed")
    return applied
