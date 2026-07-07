#!/usr/bin/env python3
"""validate_env.py - self-check kit for rbfe_map.py's OpenFE backend.

Run INSIDE your conda env (keep rbfe_map.py + test_ligands.sdf alongside):
    conda activate <your-openfe-env>
    python validate_env.py --sdf test_ligands.sdf
Paste the ENTIRE output back to the assistant.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import argparse, importlib, sys, traceback
from pathlib import Path
PASS, FAIL, WARN = "PASS", "FAIL", "WARN"
results = []
def report(stage, status, detail=""):
    results.append((stage, status))
    print(f"[{status:4}] {stage}" + (f"  -- {detail}" if detail else ""))
def stage_imports():
    print("\n=== STAGE 1: imports ===")
    ok = True
    for p in ["rdkit", "networkx", "openfe", "lomap", "konnektor", "gufe"]:
        try:
            m = importlib.import_module(p); report(f"import {p}", PASS, getattr(m, "__version__", "?"))
        except Exception as e:
            report(f"import {p}", FAIL, f"{type(e).__name__}: {e}")
            if p in ("openfe", "lomap", "gufe", "konnektor"):
                ok = False
    return ok
def stage_api():
    print("\n=== STAGE 2: API surface ===")
    ok = True
    for modname, attrs in [
        ("openfe.setup.ligand_network_planning",
         ["generate_minimal_redundant_network", "generate_minimal_spanning_network", "generate_radial_network"]),
        ("openfe.setup.atom_mapping", ["lomap_scorers"]),
        ("lomap", ["LomapAtomMapper"]), ("gufe", ["SmallMoleculeComponent"])]:
        try:
            mod = importlib.import_module(modname)
        except Exception as e:
            report(f"module {modname}", FAIL, f"{type(e).__name__}: {e}"); ok = False; continue
        for a in attrs:
            if hasattr(mod, a):
                report(f"{modname}.{a}", PASS)
            else:
                report(f"{modname}.{a}", FAIL, "attribute missing (version mismatch?)"); ok = False
    try:
        from openfe.setup.ligand_network_planning import generate_lomap_network  # noqa
        report("optional generate_lomap_network", PASS)
    except Exception:
        report("optional generate_lomap_network", WARN, "absent (ok)")
    return ok
def stage_charge_preflight(openfe_bin, charge_method):
    print("\n=== STAGE 0: OpenFE charge-generation preflight ===")
    try:
        from amber_md.gui.openfe_common import preflight_openfe_charges
    except Exception:
        # validate_env.py may be run standalone next to rbfe_map.py without the
        # package importable; fall back to a local copy of the logic.
        try:
            sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
            from amber_md.gui.openfe_common import preflight_openfe_charges
        except Exception as e:
            report("charge preflight import", WARN,
                   f"could not import preflight_openfe_charges ({e})")
            return
    pf = preflight_openfe_charges(openfe_bin=openfe_bin,
                                  charge_method=charge_method)
    for k, v in pf["info"].items():
        print(f"        {k}: {v}")
    for w in pf["warnings"]:
        report("charge preflight", WARN, w)
    for er in pf["errors"]:
        report("charge preflight", FAIL, er)
    if pf["ok"] and not pf["warnings"]:
        report("charge preflight", PASS,
               f"{charge_method} charges available")
    elif pf["ok"]:
        report("charge preflight", PASS,
               f"{charge_method} usable (see warnings)")


def _import_rbfe_map():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    return importlib.import_module("rbfe_map")
def stage_smoke(sdf, r):
    print("\n=== STAGE 3: openfe backend smoke test ===")
    try:
        c, m = r.run_openfe_backend(r.load_molecules([sdf]), "redundant", None, 2, False)
        report("openfe backend ran", PASS, f"{len(c)} edges"); return c
    except Exception as e:
        report("openfe backend ran", FAIL, f"{type(e).__name__}: {e}"); traceback.print_exc(); return None
def stage_crosscheck(sdf, r):
    print("\n=== STAGE 4: rdkit fallback cross-check ===")
    try:
        c, m = r.run_rdkit_backend(r.load_molecules([sdf]), "redundant", None, 2, False, 300, False)
        report("rdkit backend ran", PASS, f"{len(c)} edges"); return c
    except Exception as e:
        report("rdkit backend ran", FAIL, f"{type(e).__name__}: {e}"); traceback.print_exc(); return None
def stage_compare(o_e, r_e):
    print("\n=== STAGE 5: side-by-side ===")
    es = lambda E: sorted("~".join(sorted((e["lig_a"], e["lig_b"]))) for e in E) if E else []
    o, r = es(o_e), es(r_e)
    print(f"  openfe edges ({len(o)}): {o}")
    print(f"  rdkit  edges ({len(r)}): {r}")
    if o and r:
        common = set(o) & set(r)
        print(f"  shared edges        : {len(common)}/{max(len(o), len(r))}")
        report("backends agree >=50%", PASS if len(common) >= 0.5*max(len(o), len(r)) else WARN)
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--sdf", default="test_ligands.sdf")
    ap.add_argument("--openfe-bin", default="openfe",
                    help="path to the openfe executable to validate against")
    ap.add_argument("--charge-method", default="am1bcc",
                    help="partial-charge method the RBFE plan will use")
    a = ap.parse_args()
    print("="*60); print("rbfe_map.py OpenFE-backend validation kit")
    print(f"python {sys.version.split()[0]}"); print("="*60)
    stage_charge_preflight(a.openfe_bin, a.charge_method)
    imports_ok = stage_imports()
    api_ok = stage_api() if imports_ok else False
    if not imports_ok:
        report("STAGE 2 skipped", WARN, "imports failed")
    try:
        r = _import_rbfe_map(); report("import rbfe_map.py", PASS)
    except Exception as e:
        report("import rbfe_map.py", FAIL, str(e)); return _summary()
    if not Path(a.sdf).exists():
        report("test SDF present", FAIL, f"{a.sdf} not found"); return _summary()
    report("test SDF present", PASS, a.sdf)
    o_e = stage_smoke(a.sdf, r) if (imports_ok and api_ok) else None
    r_e = stage_crosscheck(a.sdf, r)
    stage_compare(o_e, r_e)
    return _summary()
def _summary():
    print("\n" + "="*60)
    nf = sum(1 for _, s in results if s == FAIL)
    nps = sum(1 for _, s in results if s == PASS)
    nw = sum(1 for _, s in results if s == WARN)
    print(f"SUMMARY: {nps} PASS, {nw} WARN, {nf} FAIL")
    print("VERDICT: OpenFE backend looks USABLE." if nf == 0
          else "VERDICT: OpenFE backend NOT ready -- paste this output back.")
    print("="*60)
    return 1 if nf else 0
if __name__ == "__main__":
    sys.exit(main())
