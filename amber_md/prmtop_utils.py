"""Lightweight prmtop parsers (no third-party deps).

Provides:
  - read_residue_labels(path)   -> list[str]
  - guess_ligand_residues(path) -> list[str]      (non-standard residues)
  - suggest_ti_masks(path)      -> dict           (timask/scmask suggestions)
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path

_STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HID","HIE","HIP","CYX","CYM","ASH","GLH","LYN",
    "NALA","NARG","NASN","NASP","NCYS","NGLU","NGLN","NGLY","NHIS","NILE",
    "NLEU","NLYS","NMET","NPHE","NPRO","NSER","NTHR","NTRP","NTYR","NVAL",
    "CALA","CARG","CASN","CASP","CCYS","CGLU","CGLN","CGLY","CHIS","CILE",
    "CLEU","CLYS","CMET","CPHE","CPRO","CSER","CTHR","CTRP","CTYR","CVAL",
    "ACE","NME","NMA","NHE",
}
_NUCLEIC = {
    "DA","DT","DG","DC","DU","A","T","G","C","U",
    "DA3","DA5","DT3","DT5","DG3","DG5","DC3","DC5",
    "RA","RU","RG","RC","RA3","RA5","RU3","RU5","RG3","RG5","RC3","RC5",
}
_WATER = {"WAT","HOH","TIP3","T3P","SOL","TIP","TP3","TP4","TIP4","SPC"}
_IONS  = {
    "NA","NA+","CL","CL-","K","K+","MG","MG2","CA","CA2","ZN","ZN2",
    "FE","FE2","FE3","MN","MN2","CU","CU1","CU2","LI","LI+","RB","CS",
    "BR","BR-","I","I-","F","F-",
    "Na+","Cl-","K+","Mg2+","Ca2+","Zn2+",
}
_LIPIDS = {"POPC","POPE","POPG","DMPC","DPPC","DOPC","CHL1","CHOL"}

_IGNORE_UPPER = {x.upper()
                 for x in (_STANDARD_AA | _NUCLEIC | _WATER | _IONS | _LIPIDS)}


def read_residue_labels(prmtop_path) -> list[str]:
    """Return the ordered list of residue labels from a prmtop.

    Tolerates the optional %COMMENT lines between %FLAG and %FORMAT, and
    handles the standard 20a4 format used for RESIDUE_LABEL.
    """
    p = Path(prmtop_path)
    if not p.exists():
        raise FileNotFoundError(p)

    labels: list[str] = []
    in_section = False
    in_data    = False
    with p.open("rt", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("%FLAG"):
                in_section = (line.strip() == "%FLAG RESIDUE_LABEL")
                in_data = False
                continue
            if not in_section:
                continue
            if line.startswith("%COMMENT"):
                continue
            if line.startswith("%FORMAT"):
                in_data = True
                continue
            if not in_data:
                continue
            if line.startswith("%"):  # safety: any other directive
                break
            # 20a4 -- four characters per name, padded.
            chunks = [line[i:i + 4].strip() for i in range(0, len(line), 4)]
            chunks = [c for c in chunks if c]
            if not chunks:
                chunks = line.split()
            labels.extend(chunks)
    return labels


def guess_ligand_residues(prmtop_path) -> list[str]:
    """Return non-standard residue names found in the prmtop, in order of
    first appearance. Duplicates removed.
    """
    labels = read_residue_labels(prmtop_path)
    seen: list[str] = []
    for name in labels:
        norm = name.strip()
        if not norm:
            continue
        if norm.upper() in _IGNORE_UPPER:
            continue
        if norm not in seen:
            seen.append(norm)
    return seen


def suggest_ti_masks(prmtop_path) -> dict:
    """Suggest timask1/timask2/scmask1/scmask2 from a prmtop.

    Returns:
        dict with keys:
          'ligands'   : list of detected non-standard residue names
          'timask1', 'timask2', 'scmask1', 'scmask2'  : suggested strings
          'mode_hint' : 'relative' (2 ligands), 'absolute' (1), or 'unknown'
          'warning'   : free-form note if anything's off
    """
    ligs = guess_ligand_residues(prmtop_path)
    out: dict = {"ligands": ligs, "warning": ""}
    if len(ligs) == 2:
        out.update(
            timask1=f":{ligs[0]}", timask2=f":{ligs[1]}",
            scmask1=f":{ligs[0]}", scmask2=f":{ligs[1]}",
            mode_hint="relative")
    elif len(ligs) == 1:
        out.update(
            timask1=f":{ligs[0]}", timask2=":MOD",
            scmask1=f":{ligs[0]}", scmask2=":MOD",
            mode_hint="absolute",
            warning=(f"Only one non-standard residue ({ligs[0]}) found. "
                     "Looks like a single-ligand system -- use absolute "
                     "mode, or rebuild with a dual-topology prmtop for "
                     "relative FEP."))
    elif len(ligs) == 0:
        out.update(
            timask1=":LIG", timask2=":MOD",
            scmask1=":LIG", scmask2=":MOD",
            mode_hint="unknown",
            warning=("No non-standard residues detected. Either your prmtop "
                     "has no ligand or its residue name is in the 'standard' "
                     "list. Edit masks manually."))
    else:
        out.update(
            timask1=f":{ligs[0]}", timask2=f":{ligs[1]}",
            scmask1=f":{ligs[0]}", scmask2=f":{ligs[1]}",
            mode_hint="relative",
            warning=(f"Found {len(ligs)} non-standard residues "
                     f"({', '.join(ligs)}). Using the first two; verify "
                     "these are the ones you want to perturb."))
    return out


if __name__ == "__main__":
    import argparse, json, sys
    p = argparse.ArgumentParser(
        description="Suggest TI masks from a prmtop's residue labels.")
    p.add_argument("prmtop", type=Path)
    p.add_argument("--json", action="store_true")
    a = p.parse_args()
    s = suggest_ti_masks(a.prmtop)
    if a.json:
        print(json.dumps(s, indent=2))
    else:
        print(f"Detected ligand residues: {s['ligands']}")
        print(f"  mode hint: {s['mode_hint']}")
        print(f"  timask1 = {s['timask1']}")
        print(f"  timask2 = {s['timask2']}")
        print(f"  scmask1 = {s['scmask1']}")
        print(f"  scmask2 = {s['scmask2']}")
        if s["warning"]:
            print(f"  WARNING: {s['warning']}", file=sys.stderr)
