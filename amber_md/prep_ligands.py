#!/usr/bin/env python3
"""prep_ligands.py — clean a Maestro SDF for use with amber_md.batch

Usage:
    python prep_ligands.py input.sdf output_dir/

What it does:
    1. Reads multi-record SDF (handles V3000 from Maestro/RDKit)
    2. Writes one V2000 SDF per molecule to output_dir/
    3. Extracts charges from Maestro tags (i_epik_Tot_Q, etc.)
    4. Writes a charges.tsv manifest

Then point the batch driver at the OUTPUT DIRECTORY:
    python -m amber_md.batch --ligands output_dir/ ...

Requires: rdkit (in your amber-md conda env: `conda install -c conda-forge rdkit`)
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

import sys, re
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print(__doc__); sys.exit(1)
    src = Path(sys.argv[1]).resolve()
    out_dir = Path(sys.argv[2]).resolve()
    if not src.exists():
        print(f"ERROR: {src} not found"); sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        print("ERROR: RDKit not installed in this env.")
        print("  conda install -c conda-forge rdkit")
        sys.exit(1)

    print(f"Reading {src}...")
    suppl = Chem.SDMolSupplier(str(src), removeHs=False, sanitize=False)

    manifest = []
    n_ok = n_skip = 0
    for i, mol in enumerate(suppl, 1):
        if mol is None:
            print(f"  record #{i}: RDKit could not parse -- SKIPPED")
            n_skip += 1
            continue
        # Sanitize (catches valence errors etc.) but don't fail on them
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception as e:
            print(f"  record #{i}: sanitize warning: {e}")

        # Name
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name").strip()
        elif mol.HasProp("ID"):
            name = mol.GetProp("ID").strip()
        else:
            name = f"MOL{i}"
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", name)[:40] or f"MOL{i}"

        # Charge: prefer Maestro/Epik tags, else compute from formal charges
        charge = None
        for tag in ("i_epik_Tot_Q", "i_user_NetCharge", "r_user_NetCharge",
                    "TOTAL_CHARGE", "CHARGE", "NetCharge"):
            if mol.HasProp(tag):
                try:
                    charge = int(float(mol.GetProp(tag)))
                    break
                except Exception:
                    pass
        if charge is None:
            charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())

        # Write V2000
        fname = out_dir / f"lig_{i:04d}_{safe}.sdf"
        w = Chem.SDWriter(str(fname))
        w.SetForceV3000(False)
        w.write(mol)
        w.close()

        # Validate
        lines = fname.read_text().splitlines()
        try:
            n_atoms = int(lines[3][0:3].strip())
        except Exception:
            n_atoms = 0
        if n_atoms == 0:
            print(f"  record #{i}: write OK but counts line bad -- SKIPPED")
            fname.unlink()
            n_skip += 1
            continue

        manifest.append((i, safe, charge, n_atoms, fname.name))
        n_ok += 1

    # Write manifest
    mfile = out_dir / "charges.tsv"
    with open(mfile, "w") as f:
        f.write("idx\tname\tcharge\tn_atoms\tfile\n")
        for row in manifest:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"\n{'='*70}")
    print(f"SUMMARY: {n_ok} ligands written to {out_dir}/")
    print(f"         {n_skip} record(s) skipped")
    print(f"Manifest: {mfile}")
    print(f"{'='*70}\n")
    print(f"{'idx':>4}  {'name':40s} {'charge':>6} {'atoms':>5}  file")
    print("-" * 90)
    for idx, name, charge, n_atoms, fname in manifest:
        print(f"{idx:>4}  {name[:40]:40s} {charge:+6d} {n_atoms:>5}  {fname}")

    print(f"\nNext: python -m amber_md.batch --ligands {out_dir} ...")

if __name__ == "__main__":
    main()
