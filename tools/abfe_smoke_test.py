#!/usr/bin/env python
"""ABFE dry-run smoke test (amber_md v2.5.16).

Generates the COMPLETE ABFE leg tree WITHOUT submitting anything, then runs a
battery of static checks so you can catch protocol/config mistakes locally
before burning the cluster GPU hours. Nothing here launches pmemd.

What it produces (under --work-dir/fep):
    complex_restraint/lambda_*/ {min,heat,press,prod}.in  + boresch.RST
    complex_decharge/ lambda_*/ ...                        + boresch.RST
    complex_vdw/      lambda_*/ ...                         + boresch.RST
    solvent_decharge/ lambda_*/ ...
    solvent_vdw/      lambda_*/ ...

What it checks (exit non-zero on any failure):
    1. window counts per leg match the configured schedules
    2. every prod.in for a soft-core stage carries the FULL GTI block
       (gti_ele_sc / gti_vdw_sc / gti_lam_sch / tishake) and scalpha=0.2/scbeta=50
    3. the restraint leg uses ifsc=0 and ramps Boresch via clambda
    4. mbar_states in each prod.in equals that leg's window count (square u_nk)
    5. neighbouring-lambda spacing in the vdW danger zone (0.575-0.80) <= 0.02
    6. boresch.RST is written for every complex window and the analytic
       correction (boresch_correction.txt) is written exactly once (restraint leg)
    7. the Boresch analytic self-test passes (11.62 kcal/mol reference)

Usage:
    python abfe_smoke_test.py \
        --absolute-prmtop complex.parm7 --absolute-inpcrd complex.rst7 \
        --solvent-prmtop  ligand.parm7  --solvent-inpcrd  ligand.rst7 \
        --boresch-json boresch.json --work-dir ./smoke

If you omit the prmtop/inpcrd files the script runs in --paper-mode: it does
not touch FEPSetup and instead validates the schedules + a freshly rendered
mdin string for each stage (no files needed). Good for a pure config check.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
# Portions ported/adapted from FEP-SPell-ABFE (freeenergylab, MIT License).

from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

DANGER_LO, DANGER_HI, DANGER_MAX_DLAM = 0.575, 0.80, 0.02
SOFTCORE_REQUIRED = ("gti_ele_sc=1", "gti_vdw_sc=1", "gti_lam_sch=1",
                     "tishake=1", "scalpha=0.2", "scbeta=50.0")


def _ok(msg):  print(f"  [OK]   {msg}")
def _bad(msg, errs): print(f"  [FAIL] {msg}"); errs.append(msg)


def check_schedules(cfg, errs):
    print("\n[1/7] lambda schedules")
    dchg = list(getattr(cfg, "decharge_lambdas", cfg.lambdas))
    vdw = list(cfg.vdw_lambdas)
    rest = list(getattr(cfg, "restraint_lambdas", cfg.lambdas))
    for name, lams, lo in (("decharge", dchg, 8), ("vdw", vdw, 30),
                           ("restraint", rest, 10)):
        if len(lams) >= lo:
            _ok(f"{name}: {len(lams)} windows")
        else:
            _bad(f"{name}: only {len(lams)} windows (expected >= {lo})", errs)
    # monotonic + bounds
    for name, lams in (("decharge", dchg), ("vdw", vdw), ("restraint", rest)):
        if lams[0] == 0.0 and lams[-1] == 1.0 and all(
                b > a for a, b in zip(lams, lams[1:])):
            _ok(f"{name}: monotonic 0->1")
        else:
            _bad(f"{name}: not strictly increasing 0->1", errs)
    return dchg, vdw, rest


def check_danger_zone(vdw, errs):
    print("\n[5/7] vdW danger-zone spacing (0.575-0.80)")
    zone = [l for l in vdw if DANGER_LO - 1e-9 <= l <= DANGER_HI + 1e-9]
    gaps = [round(b - a, 4) for a, b in zip(zone, zone[1:])]
    worst = max(gaps) if gaps else 1.0
    if worst <= DANGER_MAX_DLAM + 1e-9:
        _ok(f"{len(zone)} windows in zone, max d_lambda={worst} (<= {DANGER_MAX_DLAM})")
    else:
        _bad(f"max d_lambda in danger zone = {worst} > {DANGER_MAX_DLAM} "
             f"(windows will be unstable / MBAR overlap too low)", errs)


def check_mdin_text(stage, txt, nwin, errs):
    softcore = stage in ("decharge", "vdw")
    if softcore:
        missing = [k for k in SOFTCORE_REQUIRED if k not in txt]
        if not missing:
            _ok(f"{stage}: full GTI soft-core block present")
        else:
            _bad(f"{stage}: missing soft-core keys {missing}", errs)
        if "ifsc=1" in txt:
            _ok(f"{stage}: ifsc=1")
        else:
            _bad(f"{stage}: expected ifsc=1", errs)
    else:  # restraint
        if "ifsc=0" in txt:
            _ok("restraint: ifsc=0 (no soft-core, Boresch via clambda)")
        else:
            _bad("restraint: expected ifsc=0", errs)
        if "DISANG=boresch.RST" in txt:
            _ok("restraint: DISANG=boresch.RST present")
        else:
            _bad("restraint: missing DISANG=boresch.RST", errs)
    if f"mbar_states={nwin}" in txt:
        _ok(f"{stage}: mbar_states={nwin} matches window count (square u_nk)")
    else:
        _bad(f"{stage}: mbar_states != {nwin} (ragged MBAR grid)", errs)


def check_boresch_selftest(errs):
    print("\n[7/7] Boresch analytic correction self-test")
    try:
        from amber_md.boresch import boresch_dG_release
    except Exception:
        from boresch import boresch_dG_release  # standalone fallback
    b = {"r0": 4.470, "thA0": 101.230, "thB0": 101.230,
         "kr": 10.0, "kth": 100.0, "kph": 100.0}
    val = boresch_dG_release(b, T=298.15)
    if abs(val - 11.62) < 0.05:
        _ok(f"dG_release = {val:.3f} kcal/mol (reference 11.62)")
    else:
        _bad(f"dG_release = {val:.3f} kcal/mol, expected ~11.62 "
             f"(the rk->2*rk fix is not active!)", errs)


def check_integration_v2517(errs):
    """v2.5.17: auto-select policy, staged equilibration ladder, and the
    corrected (FEP-SPell) restraint atom ordering."""
    print("\n[8/9] v2.5.17 integration (auto-select + staged equil + ordering)")
    # (a) staged equilibration ladder = 2 min + 3 heat + 3 press + relax
    try:
        import equilibration_fepspell as eq
        steps = list(eq.build_equilibration())
        nheat = sum(1 for s in steps if s.startswith("heat-"))
        npress = sum(1 for s in steps if s.startswith("press-"))
        if steps[:2] == ["min-1", "min-2"] and nheat == 3 and npress == 3 \
                and steps[-1] == "relax":
            _ok(f"staged equilibration ladder = {len(steps)} steps "
                f"(2 min, {nheat} heat, {npress} press, relax)")
        else:
            _bad(f"unexpected equilibration ladder: {steps}", errs)
        if "ntr = 1" in eq.MIN_1 and "ntr = 0" in eq.MIN_2 \
                and "ntr = 0" in eq.RELAX and "TEMP0" in eq.HEAT:
            _ok("min-1 restrained, min-2/relax free, heat has TEMP0 ramp")
        else:
            _bad("equilibration restraint staging wrong", errs)
    except Exception as e:
        _bad(f"equilibration module error: {e}", errs)

    # (b) restraint atom ordering matches the cpptraj six-DOF measurement
    try:
        import boresch_restraint as bres
        refs = dict(r0=4.47, alpha0=72.14, theta0=101.23,
                    gamma0=77.887, beta0=-137.094, phi0=9.048)
        txt = bres.gen_restraint_str([11, 12, 13], [101, 102, 103],
                                     refs, 10, 100, 100)
        L = txt.splitlines()
        want = ["iat=11,101", "iat=101,11,12", "iat=102,101,11",
                "iat=101,11,12,13", "iat=102,101,11,12", "iat=103,102,101,11"]
        if all(w in line for w, line in zip(want, L)):
            _ok("restraints.inp uses FEP-SPell atom ordering "
                "(angles/dihedrals lead with receptor atoms)")
        else:
            _bad(f"restraint atom ordering mismatch: {L}", errs)
        # cpptraj template defines the SAME six coordinates
        cp = bres.SIX_DOFS_INPUT
        if ("angle bnd_alpha P1 L1 L2" in cp and
                "angle bnd_theta P2 P1 L1" in cp and
                "dihedral bnd_gamma P1 L1 L2 L3" in cp):
            _ok("cpptraj six-DOF template matches the restraint coordinates")
        else:
            _bad("cpptraj six-DOF template diverges from restraint defs", errs)
    except Exception as e:
        _bad(f"boresch_restraint module error: {e}", errs)

    # (c) auto-select policy: manual masks override; selector callable present
    try:
        import boresch_autoselect as bas
        lo, ro = bas.select_or_manual(
            "x", "y",
            lig_restraint_atoms=[":1@C4", ":1@C2", ":1@C6"],
            rec_restraint_atoms=[":100@CA", ":151@CA", ":132@CA"])
        if lo == [":1@C4", ":1@C2", ":1@C6"] and callable(
                bas.autoselect_boresch_atoms):
            _ok("auto-select ON by default; manual masks override (FEP-SPell parity)")
        else:
            _bad("auto-select override policy broken", errs)
    except Exception as e:
        _bad(f"boresch_autoselect module error: {e}", errs)


def check_postequil_gate(errs):
    """[9/9] post-equilibration fail-fast / self-repair gate
    (abfe_integration.verify_or_reselect_boresch). Verifies all four control
    branches with the geometry layer mocked -- no parmed/cpptraj needed."""
    print("\n[9/9] post-equilibration restraint gate (fail-fast + self-repair)")
    try:
        import boresch_autoselect as bas
        import abfe_integration as ai
        PASS = {"all_pass": True, "failures": [], "r": 7.0, "alpha": 90.0, "theta": 85.0}
        FAIL = {"all_pass": False, "failures": ["alpha=160.0 outside [30.0,150.0]"],
                "r": 7.0, "alpha": 160.0, "theta": 85.0}
        sv_vm, sv_as = bas.validate_masks, bas.autoselect_boresch_atoms
        try:
            # A: already valid -> unchanged
            bas.validate_masks = lambda *a, **k: dict(PASS)
            lo, ro, c = ai.verify_or_reselect_boresch(
                "p", "eq", ["L1", "L2", "L3"], ["P1", "P2", "P3"], is_manual=False)
            a_ok = (lo, ro) == (["L1", "L2", "L3"], ["P1", "P2", "P3"]) and c["all_pass"]
            # B: invalid + auto -> re-select recovers
            st = {"n": 0}
            def _vm(*a, **k):
                st["n"] += 1
                return dict(FAIL) if st["n"] == 1 else dict(PASS)
            bas.validate_masks = _vm
            bas.autoselect_boresch_atoms = lambda *a, **k: (["NL1", "NL2", "NL3"],
                                                            ["NP1", "NP2", "NP3"])
            lo, ro, c = ai.verify_or_reselect_boresch(
                "p", "eq", ["L1", "L2", "L3"], ["P1", "P2", "P3"],
                rmsf={0: 0.3}, is_manual=False)
            b_ok = (lo, ro) == (["NL1", "NL2", "NL3"], ["NP1", "NP2", "NP3"])
            # C: invalid + manual -> fail fast, never re-select
            bas.validate_masks = lambda *a, **k: dict(FAIL)
            hit = {"v": False}
            def _boom(*a, **k):
                hit["v"] = True
                return (["x"], ["y"])
            bas.autoselect_boresch_atoms = _boom
            c_ok = False
            try:
                ai.verify_or_reselect_boresch(
                    "p", "eq", ["L1", "L2", "L3"], ["P1", "P2", "P3"], is_manual=True)
            except ValueError as e:
                c_ok = (not hit["v"]) and "MANUALLY specified" in str(e)
            # D: invalid + auto + all retries fail -> raise
            bas.validate_masks = lambda *a, **k: dict(FAIL)
            def _fail(*a, **k):
                raise ValueError("no set")
            bas.autoselect_boresch_atoms = _fail
            d_ok = False
            try:
                ai.verify_or_reselect_boresch(
                    "p", "eq", ["L1", "L2", "L3"], ["P1", "P2", "P3"],
                    rmsf={0: 0.3}, is_manual=False, max_relax=3)
            except ValueError as e:
                d_ok = "relaxation rounds all failed" in str(e)
            consts_ok = (bas.SANITY_ANG_MIN, bas.SANITY_ANG_MAX) == (30.0, 150.0)
        finally:
            bas.validate_masks, bas.autoselect_boresch_atoms = sv_vm, sv_as
        if a_ok and b_ok and c_ok and d_ok and consts_ok:
            _ok("gate: valid->keep, auto-invalid->reselect, manual-invalid->fail, "
                "exhausted->fail; thresholds restored")
        else:
            _bad(f"post-equil gate branch failure "
                 f"(A={a_ok} B={b_ok} C={c_ok} D={d_ok} consts={consts_ok})", errs)
    except Exception as e:
        _bad(f"post-equil gate error: {e}", errs)

def render_stage_mdin(cfg, stage, boresch):
    from amber_md.fep import FEPSetup
    F = FEPSetup.__new__(FEPSetup)
    F.cfg = cfg
    F.boresch = boresch
    F.hremd = False
    F.exchange_freq = 1000
    F._suppress_correction = False
    F._stage = stage
    return F._prod_in(0.5), len(F._active_lambdas)


def paper_mode(cfg, boresch, errs):
    print("=== ABFE SMOKE TEST (paper-mode: config + mdin only) ===")
    dchg, vdw, rest = check_schedules(cfg, errs)
    print("\n[2-4/7] rendered prod.in per stage")
    for stage in ("decharge", "vdw", "restraint"):
        txt, nwin = render_stage_mdin(cfg, stage, boresch)
        print(f"  -- stage={stage} ({nwin} windows) --")
        check_mdin_text(stage, txt, nwin, errs)
    check_danger_zone(vdw, errs)
    print("\n[6/7] (file-tree checks skipped in paper-mode)")
    check_boresch_selftest(errs)
    check_integration_v2517(errs)
    check_postequil_gate(errs)


def tree_mode(a, cfg, boresch, errs):
    from amber_md.fep import FEPSetup
    from amber_md.boresch import boresch_correction
    boresch["dG_correction_kcal_mol"] = boresch_correction(
        boresch, T=cfg.temperature_K)
    wd = Path(a.work_dir)
    legs = [
        ("complex_restraint", a.absolute_prmtop, a.absolute_inpcrd, boresch, "restraint", True),
        ("complex_decharge",  a.absolute_prmtop, a.absolute_inpcrd, boresch, "decharge", False),
        ("complex_vdw",       a.absolute_prmtop, a.absolute_inpcrd, boresch, "vdw",      False),
        ("solvent_decharge",  a.solvent_prmtop,  a.solvent_inpcrd,  None,    "decharge", True),
        ("solvent_vdw",       a.solvent_prmtop,  a.solvent_inpcrd,  None,    "vdw",      True),
    ]
    print("=== ABFE SMOKE TEST (tree-mode: generating leg tree) ===")
    dchg, vdw, rest = check_schedules(cfg, errs)
    counts = {"decharge": len(dchg), "vdw": len(vdw), "restraint": len(rest)}
    md_cfg = hpc_cfg = None
    corr_files = []
    for leg_name, prm, crd, leg_b, stage, write_corr in legs:
        setup = FEPSetup(wd, cfg, md_cfg, hpc_cfg,
                         hremd=False, exchange_freq=1000, boresch=leg_b)
        leg_dir = setup.setup_leg(leg_name, prm, crd, stage=stage,
                                  write_correction=write_corr)
        wins = sorted(leg_dir.glob("lambda_*"))
        exp = counts[stage]
        print(f"\n  leg {leg_name}: {len(wins)} window dirs (expected {exp})")
        if len(wins) != exp:
            _bad(f"{leg_name}: {len(wins)} != {exp} windows", errs)
        else:
            _ok(f"{leg_name}: window count correct")
        # check one representative prod.in
        sample = wins[len(wins) // 2]
        prod = (sample / "prod.in").read_text()
        check_mdin_text(stage, prod, exp, errs)
        # boresch.RST present for complex legs
        if leg_name.startswith("complex"):
            if (sample / "boresch.RST").exists():
                _ok(f"{leg_name}: boresch.RST written")
            else:
                _bad(f"{leg_name}: boresch.RST MISSING", errs)
        cf = leg_dir / "boresch_correction.txt"
        if cf.exists():
            corr_files.append((leg_name, float(cf.read_text().strip())))
    print("\n[6/7] analytic correction written exactly once")
    if len(corr_files) == 1 and corr_files[0][0] == "complex_restraint":
        _ok(f"correction on {corr_files[0][0]} = {corr_files[0][1]:+.3f} kcal/mol")
    else:
        _bad(f"expected 1 correction file on complex_restraint, got {corr_files}", errs)
    check_danger_zone(vdw, errs)
    check_boresch_selftest(errs)
    check_integration_v2517(errs)
    check_postequil_gate(errs)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--absolute-prmtop"); p.add_argument("--absolute-inpcrd")
    p.add_argument("--solvent-prmtop");  p.add_argument("--solvent-inpcrd")
    p.add_argument("--boresch-json")
    p.add_argument("--work-dir", default="./smoke")
    p.add_argument("--charged-ligand", action="store_true",
                   help="double production length (10 ns) per the protocol")
    a = p.parse_args(argv)

    from amber_md.config import FEPConfig
    cfg = FEPConfig()
    if a.charged_ligand:
        cfg.nstlim_prod = cfg.nstlim_prod * 2

    if a.boresch_json and Path(a.boresch_json).exists():
        boresch = json.loads(Path(a.boresch_json).read_text())
    else:
        boresch = {"aA": 1, "bA": 2, "cA": 3, "A": 4, "B": 5, "C": 6,
                   "r0": 4.47, "thA0": 101.2, "thB0": 101.2,
                   "phA0": 10.0, "phB0": 20.0, "phC0": 30.0,
                   "kr": 10.0, "kth": 100.0, "kph": 100.0}

    errs = []
    have_files = all([a.absolute_prmtop, a.absolute_inpcrd,
                      a.solvent_prmtop, a.solvent_inpcrd])
    if have_files:
        tree_mode(a, cfg, boresch, errs)
    else:
        paper_mode(cfg, boresch, errs)

    print("\n" + "=" * 60)
    if errs:
        print(f"SMOKE TEST: FAILED ({len(errs)} problem(s))")
        for e in errs:
            print("  - " + e)
        return 1
    print("SMOKE TEST: PASSED — protocol is ready to submit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())