"""FEP driver entrypoint - v2.4.19.

v2.4.18:
  * NEW --resume: resubmit only unfinished windows (HPC
    maintenance/preemption recovery) and re-fire the
    analyzer + cycle-closer so ABFE_RESULT.txt completes.
    Topology rebuild is skipped (existing prmtops reused).

v2.4.14:
  * ABFE mode now ALSO submits a cycle-closer LSF job that waits on both
    per-leg analyzers (-w "done(ANA_ABS) && done(ANA_SOL)") and writes
    fep/ABFE_RESULT.txt + fep/ABFE_RESULT.json. So the pipeline completes
    without --analyze and without any Python process being alive after
    submit. This is what the GUI relies on.

v2.4.12 PATCH C (kept):
  * --mode abfe auto-builds the solvent leg, reusing GAFF params.
  * --no-solvent opts out.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from .config import HPCConfig, MDConfig, FEPConfig, SystemConfig
from .fep import (FEPSetup, FEPAnalyzer,
                  relative_binding_dG, absolute_binding_dG_from_legs)
from .logger import get_logger
from .utils import detect_ligand_resname
import csv as _csv
log = get_logger()




def _assert_mask_matches(a):
    """final46: abort if the ligand residue (a.ligand_resname) is absent from the
    built complex/solvent prmtops, i.e. timask1/scmask1 would match 0 atoms.
    Returns 0 if OK, non-zero to abort run_fep before any window is submitted."""
    rn = str(getattr(a, "ligand_resname", "") or "").strip()
    if not rn:
        log.error("ABFE mask guard: ligand_resname is empty; cannot validate "
                  "alchemical masks. Pass --ligand-resname.")
        return 4
    targets = []
    if getattr(a, "absolute_prmtop", None):
        targets.append(("complex", Path(a.absolute_prmtop)))
    if getattr(a, "solvent_prmtop", None):
        targets.append(("solvent", Path(a.solvent_prmtop)))
    if not targets:
        return 0  # nothing built (e.g. prebuilt-prmtop path handled elsewhere)
    for leg, prm in targets:
        n = _count_residue_atoms(prm, rn)
        if n is None:
            log.warning("ABFE mask guard: could not parse %s to validate "
                        "resname '%s'; skipping (run may still fail).", prm, rn)
            continue
        if n == 0:
            log.error(
                "ABFE mask guard ABORT: residue ':%s' matches 0 atoms in the "
                "%s topology (%s). The alchemical masks (timask1/scmask1/"
                "crgmask=':%s') would select NOTHING, so pmemd would treat the "
                "WHOLE system as softcore and every window would blow up. This "
                "usually means the prmtop ligand residue has a DIFFERENT name "
                "(check %%FLAG RESIDUE_LABEL). Fix: relaunch with "
                "--ligand-resname matching the prmtop, or rebuild the topology.",
                rn, leg, prm, rn)
            return 5
        log.info("ABFE mask guard OK: ':%s' matches %d atoms in %s topology.",
                 rn, n, leg)
    return 0


def _count_residue_atoms(prmtop_path, resname):
    """Return the number of atoms in residues named `resname` in a prmtop.
    Tries parmed; falls back to raw RESIDUE_LABEL/POINTERS parsing. None on
    failure."""
    try:
        import parmed
        p = parmed.load_file(str(prmtop_path))
        return sum(len(r.atoms) for r in p.residues if r.name.strip() == resname)
    except Exception:
        pass
    # Fallback: parse RESIDUE_LABEL + RESIDUE_POINTER + atom count from prmtop.
    try:
        txt = Path(prmtop_path).read_text(errors="replace")
        def _flag(name):
            i = txt.find("%FLAG " + name)
            if i < 0:
                return None
            j = txt.find("%FLAG", i + 1)
            block = txt[i:(j if j > 0 else len(txt))]
            k = block.find("%FORMAT")
            k = block.find("\n", k) + 1
            return block[k:]
        labels_b = _flag("RESIDUE_LABEL"); ptr_b = _flag("RESIDUE_POINTER")
        ptrs_a = _flag("POINTERS")
        if labels_b is None or ptr_b is None or ptrs_a is None:
            return None
        labels = labels_b.split()
        ptrs = [int(x) for x in ptr_b.split()]
        natom = int(ptrs_a.split()[0])
        total = 0
        for idx, lab in enumerate(labels):
            if lab.strip() != resname:
                continue
            start = ptrs[idx]
            end = ptrs[idx + 1] if idx + 1 < len(ptrs) else natom + 1
            total += (end - start)
        return total
    except Exception:
        return None


def _load_persisted_resname(a):
    """final46: read ligand_resname from a prior build's
    abfe_topology_inputs.json (or solvent_topology_inputs.json) under the work
    dir, so --resume reuses the SAME residue name the prmtop was built with.
    Returns the resname string or None if not found."""
    try:
        wd = Path(a.work_dir).expanduser()
    except Exception:
        return None
    for rel in ("abfe_topology_inputs.json",
                "build/abfe_topology_inputs.json",
                "solvent_topology_inputs.json",
                "build_solvent/solvent_topology_inputs.json"):
        p = wd / rel
        if p.exists():
            try:
                rn = json.loads(p.read_text()).get("ligand_resname")
                if rn:
                    return str(rn)
            except Exception:
                continue
    return None


def _derive_abfe_masks(a):
    """For ABFE every alchemical mask must reference the ligand residue. Derive
    from a.ligand_resname UNLESS the user explicitly set a mask. Single-topology
    decoupling: region1 = :<resname>, region2 = '' (empty)."""
    lig = ":" + str(a.ligand_resname).lstrip(":")
    DEF = {"timask1": ":LIG", "timask2": ":MOD",
           "scmask1": ":LIG", "scmask2": ":MOD"}
    if getattr(a, "timask1", DEF["timask1"]) == DEF["timask1"]: a.timask1 = lig
    if getattr(a, "timask2", DEF["timask2"]) == DEF["timask2"]: a.timask2 = ""
    if getattr(a, "scmask1", DEF["scmask1"]) == DEF["scmask1"]: a.scmask1 = lig
    if getattr(a, "scmask2", DEF["scmask2"]) == DEF["scmask2"]: a.scmask2 = ""
    if not getattr(a, "crgmask", None): a.crgmask = lig
    return a


def _build_configs(a):
    md  = MDConfig(temperature_K=a.temperature)
    # PATCH Bug 6 (v2.4.19): thread per-stage lambda schedules through to
    # FEPConfig. Previously --lambdas was silently ignored for ABFE legs,
    # which used the FEPConfig stage-schedule defaults.
    fep_kwargs = dict(
        lambdas=a.lambdas,
        nstlim_eq=a.nstlim_eq, nstlim_prod=a.nstlim_prod,
        restraint_nstlim_prod=getattr(a, 'restraint_nstlim_prod', 2_000_000),
        dt_ps=a.dt, cutoff_A=a.cutoff,
        ntpr=a.ntpr, ntwx=a.ntwx,
        temperature_K=a.temperature,
        timask1=a.timask1, timask2=a.timask2,
        scmask1=a.scmask1, scmask2=a.scmask2,
    )
    if getattr(a, "crgmask", None):
        fep_kwargs["crgmask"] = a.crgmask
    if getattr(a, "decharge_lambdas", None):
        fep_kwargs["decharge_lambdas"] = tuple(a.decharge_lambdas)
    if getattr(a, "vdw_lambdas", None):
        fep_kwargs["vdw_lambdas"] = tuple(a.vdw_lambdas)
    if getattr(a, "fine_restraint", False):
        fep_kwargs["use_fine_restraint_lambdas"] = True
    fep = FEPConfig(**fep_kwargs)
    # v2.5.62: hard guard -- refuse to launch with a >1 fs production dt.
    from ._dt_guard import assert_prod_dt_safe
    assert_prod_dt_safe(fep)
    assert_prod_dt_safe(md)
    # Default project comes from HPCConfig (your-project). Only override
    # it when the user explicitly passed --project, so FEP jobs never get
    # submitted under the bogus "-P default" project.
    hpc_kwargs = dict(
        queue_gpu=a.queue, walltime=a.walltime, n_gpu=a.n_gpu,
        modules=a.modules, venv_activate=a.venv,
    )
    if getattr(a, "project", None):
        hpc_kwargs["project"] = a.project
    hpc = HPCConfig(**hpc_kwargs)
    if getattr(a, "queue_cpu", None):
        hpc.queue_cpu = a.queue_cpu
    return md, fep, hpc


def _maybe_build_abfe_topology(a, wd):
    a._complex_pdb = None; a._ligand_mol2 = None; a._ligand_frcmod = None
    if a.mode != "abfe": return 0
    if a.absolute_prmtop and a.absolute_inpcrd: return 0
    if not (a.protein_pdb and a.ligand_file):
        log.error("--mode abfe requires --protein-pdb AND --ligand-file."); return 1
    from .abfe_topology import build_abfe_topology
    try:
        prm, crd, cpx_pdb = build_abfe_topology(
            protein_pdb=a.protein_pdb,
            ligand_file=a.ligand_file,
            ligand_index=a.ligand_index,
            ligand_resname=a.ligand_resname,
            work_dir=wd,
            sys_cfg=SystemConfig(ligand_charge=a.ligand_charge,
                                 charge_method=a.charge_method,
                                 box_buffer_A=a.box_buffer),
        )
        a.absolute_prmtop = prm; a.absolute_inpcrd = crd; a._complex_pdb = cpx_pdb
        ij = wd / "abfe_topology_inputs.json"
        if ij.exists():
            d = json.loads(ij.read_text())
            a._ligand_mol2 = d.get("ligand_mol2")
            a._ligand_frcmod = d.get("ligand_frcmod")
        return 0
    except Exception as e:
        log.exception("ABFE topology build failed: %s", e); return 2


def _maybe_build_solvent_topology(a, wd):
    if a.mode != "abfe": return 0
    if a.solvent_prmtop and a.solvent_inpcrd: return 0
    if not a.ligand_file and not (a._ligand_mol2 and a._ligand_frcmod):
        log.error("Cannot build solvent leg: no ligand info."); return 1
    from .abfe_topology import build_solvent_only_topology
    try:
        s_prm, s_crd = build_solvent_only_topology(
            ligand_file=a.ligand_file,
            ligand_index=a.ligand_index,
            ligand_resname=a.ligand_resname,
            work_dir=wd,
            sys_cfg=SystemConfig(ligand_charge=a.ligand_charge,
                                 charge_method=a.charge_method,
                                 box_buffer_A=a.box_buffer),
            reuse_mol2=Path(a._ligand_mol2) if a._ligand_mol2 else None,
            reuse_frcmod=Path(a._ligand_frcmod) if a._ligand_frcmod else None,
        )
        a.solvent_prmtop = s_prm; a.solvent_inpcrd = s_crd
        return 0
    except Exception as e:
        log.exception("ABFE solvent build failed: %s", e); return 2


def _read_edges_csv(path):
    rows = []
    with open(path, newline="") as f:
        for row in _csv.DictReader(f):
            rows.append(row)
    return rows


def _extract_named_record(ligands_file, name, dest_dir):
    """Pull a single named ligand record from a multi-record SDF/MOL2 into its
    own single-record file (LigandParametrizer is single-record only)."""
    from pathlib import Path as _P
    dest_dir = _P(dest_dir); dest_dir.mkdir(parents=True, exist_ok=True)
    s = _P(ligands_file); text = s.read_text(errors="replace"); suf = s.suffix.lower()
    if suf in (".sdf", ".mol"):
        for rec in text.split("$$$$"):
            rec_s = rec.strip("\n")
            if not rec_s.strip():
                continue
            if rec_s.splitlines()[0].strip().lower() == str(name).lower():
                out = dest_dir / f"{name}.sdf"; out.write_text(rec_s + "\n$$$$\n")
                return out
        raise ValueError(f"Ligand '{name}' not found in {s}")
    if suf == ".mol2":
        for blk in text.split("@<TRIPOS>MOLECULE"):
            if not blk.strip():
                continue
            if blk.strip().splitlines()[0].strip().lower() == str(name).lower():
                out = dest_dir / f"{name}.mol2"
                out.write_text("@<TRIPOS>MOLECULE" + blk)
                return out
        raise ValueError(f"Ligand '{name}' not found in {s}")
    raise ValueError(f"Unsupported ligand container format: {suf}")


def _run_single_edge(a, edge_name, lig_a, lig_b, edge_wd):
    from .rbfe_topology import build_rbfe_edge_topology
    from .config import SystemConfig
    import copy
    edge_wd = Path(edge_wd).expanduser().resolve(); edge_wd.mkdir(parents=True, exist_ok=True)
    recs = edge_wd / "ligand_records"
    a_file = _extract_named_record(a.ligands, lig_a, recs)
    b_file = _extract_named_record(a.ligands, lig_b, recs)
    log.info("=== RBFE edge %s : building dual-topology ===", edge_name)
    topo = build_rbfe_edge_topology(
        protein_pdb=a.protein_pdb, ligand_a_file=a_file, ligand_b_file=b_file,
        work_dir=edge_wd / "topo",
        sys_cfg=SystemConfig(ligand_charge=a.ligand_charge,
                             charge_method=a.charge_method,
                             box_buffer_A=a.box_buffer),
        auto_protonation=True)
    ea = copy.copy(a)
    ea.mode = "legacy"; ea.work_dir = edge_wd
    ea.complex_prmtop = Path(topo["complex_prmtop"]); ea.complex_inpcrd = Path(topo["complex_inpcrd"])
    ea.solvent_prmtop = Path(topo["solvent_prmtop"]); ea.solvent_inpcrd = Path(topo["solvent_inpcrd"])
    ea.absolute_prmtop = None; ea.absolute_inpcrd = None
    ea.auto_boresch = False; ea.boresch_json = None
    ea.timask1 = topo["timask1"]; ea.timask2 = topo["timask2"]
    ea.scmask1 = topo["scmask1"]; ea.scmask2 = topo["scmask2"]
    ea.ligand_resname = "L1"
    return run_fep(ea)


def run_rbfe(a):
    """End-to-end Amber RBFE campaign from an rbfe_map network."""
    wd = Path(a.work_dir).expanduser().resolve(); wd.mkdir(parents=True, exist_ok=True)
    if not a.ligands:
        log.error("--mode rbfe requires --ligands."); return 1
    if not a.protein_pdb:
        log.error("--mode rbfe requires --protein-pdb."); return 1
    edges_csv = a.edges_csv
    if edges_csv is None:
        map_dir = a.rbfe_map_dir or (wd / "rbfe_map")
        edges_csv = Path(map_dir) / "edges.csv"
    edges_csv = Path(edges_csv)
    if not edges_csv.exists():
        log.error("RBFE edges.csv not found: %s. Run rbfe_map first (or pass "
                  "--edges-csv).", edges_csv); return 2
    rows = _read_edges_csv(edges_csv)
    if not rows:
        log.error("edges.csv is empty: %s", edges_csv); return 2
    want = set(a.only_edges) if a.only_edges else None
    edges_dir = wd / "edges"; edges_dir.mkdir(parents=True, exist_ok=True)
    summary = {"edges": {}, "edges_csv": str(edges_csv)}
    n_ok = n_fail = 0
    for row in rows:
        edge = row.get("edge") or f"{row.get('lig_a')}~{row.get('lig_b')}"
        if want is not None and edge not in want:
            continue
        lig_a = row.get("lig_a"); lig_b = row.get("lig_b")
        if not lig_a or not lig_b:
            log.warning("RBFE edge %s missing lig_a/lig_b; skipping.", edge); continue
        edge_wd = edges_dir / edge
        try:
            rc = _run_single_edge(a, edge, lig_a, lig_b, edge_wd)
            summary["edges"][edge] = {"rc": rc, "work_dir": str(edge_wd)}
            n_ok += (rc == 0); n_fail += (rc != 0)
            if rc != 0:
                log.error("RBFE edge %s returned rc=%s", edge, rc)
        except Exception as e:  # noqa: BLE001
            n_fail += 1
            summary["edges"][edge] = {"rc": "exception", "error": str(e),
                                      "work_dir": str(edge_wd)}
            log.exception("RBFE edge %s failed: %s", edge, e)
    (wd / "RBFE_CAMPAIGN.json").write_text(json.dumps(summary, indent=2))
    log.info("RBFE campaign: %d edge(s) ok, %d failed. Per-edge ddG in each "
             "edge's fep/summary.json once windows finish.", n_ok, n_fail)
    log.info("RBFE campaign manifest: %s", wd / "RBFE_CAMPAIGN.json")
    return 0 if n_fail == 0 else 3


def run_fep(a):
    # Resolve the ligand residue name BEFORE building configs, because the ABFE
    # alchemical masks are derived from it. Priority: explicit --ligand-resname
    # > auto-detect from .mol2 > "LIG".
    if getattr(a, "mode", "legacy") == "abfe":
        if not getattr(a, "ligand_resname", None):
            # final46 FIX (resume-safe resname): the alchemical masks MUST match
            # the residue name baked into the already-built prmtop. On a --resume
            # sweep the user typically passes neither --ligand-resname nor
            # --ligand-file, which previously fell through to the "LIG" default
            # and silently rewrote :UNK masks to :LIG -> 0-atom masks -> the
            # whole system became softcore and every window blew up. Reload the
            # persisted resname from abfe_topology_inputs.json FIRST.
            _persisted_rn = _load_persisted_resname(a)
            if _persisted_rn:
                a.ligand_resname = _persisted_rn
                log.info("ligand resname: %s (reused from previous build)",
                         _persisted_rn)
            elif getattr(a, "ligand_file", None):
                rn, src = detect_ligand_resname(a.ligand_file,
                                                getattr(a, "ligand_index", 0))
                a.ligand_resname = rn
                log.info("ligand resname: %s (auto-detected from %s)", rn, src)
            else:
                a.ligand_resname = "LIG"
                log.info("ligand resname: LIG (default; no ligand file)")
        else:
            log.info("ligand resname: %s (user-specified)", a.ligand_resname)
        _derive_abfe_masks(a)
        log.info("ABFE masks -> timask1='%s' timask2='%s' scmask1='%s' "
                 "scmask2='%s' crgmask='%s'", a.timask1, a.timask2,
                 a.scmask1, a.scmask2, getattr(a, "crgmask", None))
    else:
        if not getattr(a, "ligand_resname", None):
            a.ligand_resname = "LIG"
    md_cfg, fep_cfg, hpc_cfg = _build_configs(a)
    wd = Path(a.work_dir).expanduser().resolve()
    wd.mkdir(parents=True, exist_ok=True)

    # PATCH Bug 2 (v2.4.19): honour "--resume skips topology rebuild".
    # Adopt existing prmtops/boresch.json on disk so the _maybe_build_*
    # guards short-circuit and Boresch atoms/correction are NOT re-derived.
    if a.resume:
        cpx_prm = wd / "build" / "complex.prmtop"
        cpx_crd = wd / "build" / "complex.inpcrd"
        sol_prm = wd / "build_solvent" / "solvent.prmtop"
        sol_crd = wd / "build_solvent" / "solvent.inpcrd"
        bj      = wd / "boresch.json"
        if cpx_prm.exists() and cpx_crd.exists() and not a.absolute_prmtop:
            a.absolute_prmtop, a.absolute_inpcrd = cpx_prm, cpx_crd
            log.info("  resume: reusing complex topology %s", cpx_prm)
        if sol_prm.exists() and sol_crd.exists() and not a.solvent_prmtop:
            a.solvent_prmtop, a.solvent_inpcrd = sol_prm, sol_crd
            log.info("  resume: reusing solvent topology %s", sol_prm)
        if bj.exists() and not a.boresch_json:
            a.boresch_json = bj
            log.info("  resume: reusing Boresch restraints %s", bj)

    rc = _maybe_build_abfe_topology(a, wd)
    if rc != 0: return rc
    if not a.no_solvent:
        rc = _maybe_build_solvent_topology(a, wd)
        if rc != 0: return rc

    # final46 FIX (fail-fast mask guard): verify the alchemical region actually
    # selects atoms in the built prmtop BEFORE submitting 76 doomed windows.
    # The :LIG-vs-:UNK mismatch (resume default bug) made timask1 match 0 atoms,
    # which pmemd silently treats as "whole system softcore" and runs garbage.
    if a.mode == "abfe":
        rc = _assert_mask_matches(a)
        if rc != 0: return rc

    boresch = None
    if a.absolute_prmtop and a.absolute_inpcrd:
        if a.boresch_json:
            boresch = json.loads(Path(a.boresch_json).read_text())
        elif a.auto_boresch:
            if not a.auto_boresch_pdb:
                if a.mode == "abfe" and a._complex_pdb is not None:
                    a.auto_boresch_pdb = a._complex_pdb
                else:
                    guess = wd / "build" / "complex.pdb"
                    if guess.exists(): a.auto_boresch_pdb = guess
                    else:
                        log.error("--auto-boresch needs complex PDB."); return 1
            from .boresch import select_boresch_atoms
            try:
                boresch = select_boresch_atoms(
                    a.auto_boresch_pdb, a.ligand_resname)
            except ValueError as e:
                log.error("Boresch selection failed: %s", e); return 1
        else:
            log.error("Absolute leg requested but no boresch source."); return 1
        # v2.5.17 (caveat #1 fix): bridge the legacy dict-style Boresch atoms
        # (1-based global atom serials aA/bA/cA + A/B/C, as produced by
        # select_boresch_atoms and used by the GUI --auto-boresch path) to Amber
        # @serial masks, so the in-job post-equilibration gate (B) can
        # re-validate the SAME six atoms against eq.rst. Restraint geometry:
        # bond aA(rec)-A(lig); so L1=A,L2=B,L3=C (ligand) and P1=aA,P2=bA,P3=cA
        # (receptor). Purely additive: never overwrites masks already present.
        if isinstance(boresch, dict) and not boresch.get("lig_masks"):
            try:
                _need = ("A", "B", "C", "aA", "bA", "cA")
                if all(k in boresch for k in _need):
                    boresch["lig_masks"] = [f"@{int(boresch['A'])}",
                                            f"@{int(boresch['B'])}",
                                            f"@{int(boresch['C'])}"]
                    boresch["rec_masks"] = [f"@{int(boresch['aA'])}",
                                            f"@{int(boresch['bA'])}",
                                            f"@{int(boresch['cA'])}"]
                    log.info("Boresch dict->mask bridge: lig=%s rec=%s",
                             boresch["lig_masks"], boresch["rec_masks"])
            except Exception as e:
                log.warning("Boresch dict->mask bridge skipped (%s)", e)
        # v2.5.18: ensure canonical FEP-SPell six-DOF labels are present so
        # the restraint that is SIMULATED, MEASURED and CORRECTED all use the
        # same internal coordinates (also for a boresch.json loaded from disk).
        try:
            from .boresch import canonical_dofs_from_legacy
            if isinstance(boresch, dict) and "alpha0" not in boresch \
                    and all(k in boresch for k in ("A","B","C","aA","bA","cA")):
                boresch.update(canonical_dofs_from_legacy(boresch))
        except Exception as e:
            log.warning("Boresch canonical-DOF mapping skipped (%s)", e)
        # v2.5.18: analytic standard-state correction via Deng & Roux Eq.38/40
        # (FEP-SPell analytic.py), using the SAME force constants written to
        # boresch.RST. Replaces the prior closed form that was evaluated on a
        # geometry inconsistent with the simulated restraint.
        from .boresch import boresch_correction_dengroux
        boresch["dG_correction_kcal_mol"] = boresch_correction_dengroux(
            boresch, T=a.temperature)
        log.info("Boresch correction (Deng-Roux) = %+0.3f kcal/mol",
                 boresch["dG_correction_kcal_mol"])
        (wd / "boresch.json").write_text(json.dumps(boresch, indent=2))

        # v2.5.17 (A): SETUP-TIME geometric pre-check of the Boresch atoms.
        try:
            from .boresch import _precheck_boresch_dict
            ok, why = _precheck_boresch_dict(boresch)
            if not ok:
                log.error("Boresch atoms failed setup-time geometry check: %s.", why)
                return 1
            log.info("Boresch setup-time geometry check: PASS")
        except Exception as e:
            log.warning("Boresch setup-time pre-check skipped (%s)", e)

    two_stage = getattr(fep_cfg, "two_stage", True) and a.mode == "abfe"
    # PATCH Bug 6 (v2.4.19): make the actual ABFE schedule explicit in logs.
    if two_stage:
        log.info("ABFE two-stage schedules: decharge=%d windows, vdw=%d "
                 "windows (--lambdas=%d applies to RBFE/non-two-stage legs "
                 "only).",
                 len(getattr(fep_cfg, "decharge_lambdas", fep_cfg.lambdas)),
                 len(getattr(fep_cfg, "vdw_lambdas", fep_cfg.lambdas)),
                 len(fep_cfg.lambdas))
    # tuple: (leg_name, prm, crd, boresch, stage, write_correction)
    #
    # v2.5.37 (Option A): the restraint-removal leg runs on the REAL complex
    # topology (a.absolute_prmtop) with a SINGLE ligand copy and a lambda-scaled
    # Boresch potential (k=lambda*k_full, see fep._write_boresch_RST). The old
    # dual-copy topology was removed: its :2 'dummy' copy was parameter-identical
    # to :1 (not decoupled), so the :1->:2 TI morph produced a 1/r^12 end-state
    # singularity (dV/dl ~ -1200). Original Boresch indices already point at the
    # real ligand/protein, so no remap is needed.
    restraint_prm = restraint_crd = None
    restraint_boresch = boresch
    if (two_stage and a.absolute_prmtop and a.absolute_inpcrd
            and getattr(fep_cfg, "build_restraint_topology", False)):
        log.error("build_restraint_topology=True is no longer supported "
                  "(v2.5.37 removed the dual-copy module; its :2 copy was not "
                  "decoupled). Using single-copy Option A restraint leg instead.")
    legs = []
    if a.mode == "legacy" or not two_stage:
        if a.complex_prmtop and a.complex_inpcrd:
            legs.append(("complex",  a.complex_prmtop,  a.complex_inpcrd,  None, None, True))
        if a.absolute_prmtop and a.absolute_inpcrd:
            legs.append(("absolute", a.absolute_prmtop, a.absolute_inpcrd, boresch, None, True))
        if a.solvent_prmtop and a.solvent_inpcrd:
            legs.append(("solvent",  a.solvent_prmtop,  a.solvent_inpcrd,  None, None, True))
    else:
        if a.absolute_prmtop and a.absolute_inpcrd:
            # v2.5.16: formal ABFE cycle. decharge + vdw carry the Boresch
            # POTENTIAL (so the ligand stays in the pocket while decoupling),
            # and a DEDICATED restraint-removal leg turns the Boresch
            # restraint off on the fully-interacting complex. The analytic
            # standard-state correction (+dG) is written ONCE, on the
            # restraint leg, and the cycle-closer adds it to the MD restraint
            # free energy. (Previously the correction was written on
            # complex_vdw with no MD restraint leg -- formally incomplete.)
            # v2.5.23 (Option B): restraint leg uses the dual-copy topology and
            # remapped Boresch atoms when available; falls back to complex.prmtop.
            _rprm = restraint_prm or a.absolute_prmtop
            _rcrd = restraint_crd or a.absolute_inpcrd
            # v2.5.68: the restraint MD leg is OFF by default. Its sampled FE is
            # ~0 (removing a well-centered Boresch restraint on the bound,
            # fully-interacting complex), so we skip the leg and write the
            # analytic Boresch standard-state term onto complex_vdw instead
            # (Boresch 2003 / FEP+ / OpenFE convention). Use --restraint-leg to
            # restore the explicit MD leg as a BAT.py-style validation.
            # CRITICAL INVARIANT: the Boresch POTENTIAL stays HELD ON during
            # decharge+vdw in BOTH modes (fep.py: fixed full-k restraint), so
            # the decoupled ligand cannot drift regardless of this flag.
            _want_restraint_leg = getattr(a, "restraint_leg", False)
            if _want_restraint_leg:
                legs.append(("complex_restraint", _rprm, _rcrd, restraint_boresch, "restraint", True))
                _vdw_writes_corr = False
            else:
                _vdw_writes_corr = True   # analytic Boresch term lands on complex_vdw
            legs.append(("complex_decharge",  a.absolute_prmtop, a.absolute_inpcrd, boresch, "decharge", False))
            legs.append(("complex_vdw",       a.absolute_prmtop, a.absolute_inpcrd, boresch, "vdw",      _vdw_writes_corr))
        if a.solvent_prmtop and a.solvent_inpcrd:
            legs.append(("solvent_decharge", a.solvent_prmtop, a.solvent_inpcrd, None, "decharge", True))
            legs.append(("solvent_vdw",      a.solvent_prmtop, a.solvent_inpcrd, None, "vdw",      True))
    if not legs:
        log.error("No legs to set up."); return 1

    closer_setup = FEPSetup(
        wd, fep_cfg, md_cfg, hpc_cfg,
        hremd=a.hremd, exchange_freq=a.exchange_freq, boresch=None)

    job_ids = {}
    analyzer_jids = {}
    for leg_name, prm, crd, leg_boresch, stage, write_corr in legs:
        log.info("=== Setting up leg: %s (stage=%s, corr=%s) ===",
                 leg_name, stage, write_corr)
        setup = FEPSetup(
            wd, fep_cfg, md_cfg, hpc_cfg,
            hremd=a.hremd, exchange_freq=a.exchange_freq,
            boresch=leg_boresch)
        leg_dir = setup.setup_leg(leg_name, prm, crd, stage=stage,
                                  write_correction=write_corr)
        log.info("  wrote %d lambda windows under %s",
                 len(setup._active_lambdas), leg_dir)  # PATCH: per-stage count
        # v2.5.43: explicit pre-submission Option A guard for the restraint leg.
        # setup_leg() already raises on a TI-keyword leak; here we ALSO scan
        # explicitly so batch (RBFE) runs get a clear per-edge log line and a
        # graceful abort (return 2) instead of an uncaught exception. Covers
        # every edge because RBFE routes each edge through run_fep().
        if stage == "restraint":
            try:
                from amber_md.fep import FEPSetup as _FS
                _bad = []
                for _mdin in sorted(leg_dir.glob("lambda_*/*.in")):
                    _hits = [k for k in _FS._RESTRAINT_FORBIDDEN_KW
                             if k in _mdin.read_text()]
                    if _hits:
                        _bad.append("%s: %s" %
                                    (_mdin.relative_to(leg_dir), ", ".join(_hits)))
                if _bad:
                    log.error("restraint-input guard: TI keyword(s) leaked into "
                              "Option A leg %s -- NOT submitting:\n  %s",
                              leg_name, "\n  ".join(_bad))
                    log.error("This is a generator bug (a stage template still "
                              "emits TI keywords). Aborting run_fep.")
                    return 2
                log.info("  restraint-input guard: OK (Option A, no TI keywords) "
                         "for leg %s", leg_name)
            except Exception as _rg:
                log.warning("restraint-input guard could not run (%s); "
                            "continuing.", _rg)
        # v2.5.31b: TWO-PHASE GPU-free mdin gate (content pre-script, then -ref
        # post-script) -- both before any bsub. See mdin_validator.validate_leg.
        _do_check = not getattr(a, "skip_mdin_check", False)
        if _do_check:
            try:
                from amber_md.mdin_validator import validate_leg, format_issues
                _e, _w = validate_leg(leg_dir, require_script=False)
                if _w:
                    log.warning("mdin validator (content): %d warning(s) in %s:\n%s",
                                len(_w), leg_name, format_issues([], _w))
                if _e:
                    log.error("mdin validator (content): %d FATAL issue(s) in %s -- "
                              "NOT submitting:\n%s", len(_e), leg_name,
                              format_issues(_e, []))
                    log.error("Fix the generator or rerun with --skip-mdin-check "
                              "(not recommended). Aborting run_fep.")
                    return 2
                log.info("  mdin validator (content): OK for leg %s", leg_name)
            except Exception as _ve:
                log.warning("mdin content validator could not run (%s); continuing.",
                            _ve)

        # Write run + analyze scripts now so the -ref check has a real file.
        if a.hremd:
            setup.build_lsf_hremd(leg_dir, leg_name)
        else:
            setup.build_lsf_array(leg_dir, leg_name)
        setup._build_analyze_lsf(leg_dir, leg_name)

        if _do_check:
            try:
                from amber_md.mdin_validator import validate_leg, format_issues
                _e2, _w2 = validate_leg(leg_dir, require_script=True)
                _se = [x for x in _e2 if x.where.endswith(".lsf") or "run_" in x.where]
                _sw = [x for x in _w2 if x.where.endswith(".lsf")
                       or "run_" in x.where or "no run_" in x.msg]
                if _sw:
                    log.warning("mdin validator (-ref): %d warning(s) in %s:\n%s",
                                len(_sw), leg_name, format_issues([], _sw))
                if _se:
                    log.error("mdin validator (-ref): %d FATAL run-script issue(s) in "
                              "%s -- NOT submitting:\n%s", len(_se), leg_name,
                              format_issues(_se, []))
                    log.error("Aborting run_fep (use --skip-mdin-check to override).")
                    return 2
                log.info("  mdin validator (-ref): OK for leg %s", leg_name)
            except Exception as _ve:
                log.warning("mdin -ref validator could not run (%s); continuing.", _ve)

        if a.submit or a.resume:
            if a.resume:
                jid = setup.submit_leg_resume(leg_dir, leg_name)
            else:
                jid = setup.submit_leg(leg_dir, leg_name)
            job_ids[leg_name] = jid
            analyzer_jids[leg_name] = jid.get("analyze")
            log.info("  %s MD=%s analyze=%s",
                     "resumed" if a.resume else "submitted",
                     jid.get("md"), jid.get("analyze"))

    # v2.4.14: cycle-closer for ABFE
    if a.mode == "abfe":
        cc_script = closer_setup.build_cycle_closer_lsf(wd / "fep")
        log.info("Wrote cycle-closer script: %s", cc_script)
        need = (["complex_restraint","complex_decharge","complex_vdw",
                 "solvent_decharge","solvent_vdw"]
                if getattr(fep_cfg, "two_stage", True)
                else ["absolute","solvent"])
        if (a.submit or a.resume) and all(k in analyzer_jids for k in need):
            cc_jid = closer_setup.submit_cycle_closer(
                wd / "fep", analyzer_jids)
            if cc_jid:
                job_ids["cycle_closer"] = {"md": None, "analyze": cc_jid}
                log.info("ABFE_RESULT.txt will appear in %s when both legs "
                         "finish.", wd / "fep")

    if job_ids:
        (wd / "fep" / "job_ids.json").write_text(json.dumps(job_ids, indent=2))

    if a.analyze:
        # v2.5.15: GATE -> RECOVER (logged) -> ANALYZE. Before analyzing any
        # leg we ensure every window has a complete prod.out, attempting
        # recovery for any that do not, and we record every recovery action to
        # fep/<leg>/recovery_log.json + fep/RECOVERY_REPORT.txt so the steps
        # taken can be troubleshot later. Disable with --no-auto-recover.
        results = {}
        _auto_recover = (a.mode == "abfe"
                         and not getattr(a, "no_auto_recover", False))
        _leg_records = []
        _rerun = None
        if _auto_recover:
            try:
                from .recovery import ensure_leg_complete, write_recovery_report
                from .abfe_self_heal_cli import (
                    make_local_rerun, make_bsub_rerun, make_dry_rerun)
                if a.self_heal_mode == "local":
                    _rerun = make_local_rerun()
                elif a.self_heal_mode == "dry-run":
                    _rerun = make_dry_rerun()
                else:
                    _rerun = make_bsub_rerun(
                        a.project or getattr(fep_cfg, "project",
                                             "your-project"),
                        a.queue, a.walltime,
                        getattr(fep_cfg, "fep_mem_mb", 8192))
            except Exception as e:  # noqa: BLE001
                log.error("auto-recovery unavailable (%s); analyzing as-is", e)
                _auto_recover = False
        for leg_name, _, _, _, _stage, _ in legs:  # PATCH: legs are 6-tuples
            leg_dir = wd / "fep" / leg_name
            # per-leg lambda schedule (decharge/vdw/common)
            if _stage == "decharge":
                _lams = (a.decharge_lambdas
                         or getattr(fep_cfg, "decharge_lambdas", fep_cfg.lambdas))
            elif _stage == "vdw":
                _lams = (a.vdw_lambdas
                         or getattr(fep_cfg, "vdw_lambdas", fep_cfg.lambdas))
            else:
                _lams = fep_cfg.lambdas
            if _auto_recover and leg_dir.exists():
                log.info("=== recovery gate: %s ===", leg_name)
                _rec = ensure_leg_complete(leg_dir, _lams, _rerun,
                                           getattr(a, "self_heal_attempts", 3),
                                           logger=log)
                _leg_records.append(_rec)
                if not _rec.get("complete", True):
                    log.error("  %s NOT fully recovered; analyzing for "
                              "diagnostics only (result NOT trusted).",
                              leg_name)
            res = FEPAnalyzer(leg_dir, fep_cfg.lambdas,
                              temperature_K=a.temperature).run()
            results[leg_name] = res
            (leg_dir / "summary.json").write_text(
                json.dumps(res, indent=2, default=str))
        if "absolute" in results and "solvent" in results:
            dG_bind = absolute_binding_dG_from_legs(
                results["absolute"], results["solvent"])
            if dG_bind is not None:
                dG_c = results["absolute"]["dG_kcal_mol"]
                dG_s = results["solvent"]["dG_kcal_mol"]
                msg = ("\n================ ABFE RESULT ================\n"
                       f"  dG_complex+restr  = {dG_c:+8.3f} kcal/mol\n"
                       f"  dG_solvent        = {dG_s:+8.3f} kcal/mol\n"
                       f"  dG_bind           = {dG_bind:+8.3f} kcal/mol\n"
                       "=============================================\n")
                log.info(msg)
                (wd / "fep" / "ABFE_RESULT.txt").write_text(msg)
        (wd / "fep" / "summary.json").write_text(
            json.dumps(results, indent=2, default=str))
        if _auto_recover and _leg_records:
            try:
                from .recovery import write_recovery_report
                write_recovery_report(wd / "fep", _leg_records,
                                      logger=log)
            except Exception as e:  # noqa: BLE001
                log.error("could not write recovery report: %s", e)


    # --- self-heal (final61): guarantee a complete prod.out for every window of
    #     every leg. Diagnoses missing/incomplete windows, adjusts their mdin
    #     parameters, and re-runs. Opt-in via --self-heal. ---
    if getattr(a, "self_heal", False) and a.mode == "abfe":
        try:
            from .abfe_self_heal import heal_leg
            from .abfe_self_heal_cli import (
                make_local_rerun, make_bsub_rerun, make_dry_rerun)
        except Exception as e:
            log.error("self-heal unavailable: %s", e)
        else:
            if a.self_heal_mode == "local":
                rerun = make_local_rerun()
            elif a.self_heal_mode == "dry-run":
                rerun = make_dry_rerun()
            else:
                rerun = make_bsub_rerun(
                    a.project or getattr(fep_cfg, "project",
                                         "your-project"),
                    a.queue, a.walltime,
                    getattr(fep_cfg, "fep_mem_mb", 8192))
            log.info("=== self-heal: ensuring prod.out for %d leg(s) ===",
                     len(legs))
            for leg_name, _, _, _, stage, _ in legs:
                leg_dir = wd / "fep" / leg_name
                if not leg_dir.exists():
                    continue
                if stage == "decharge":
                    lams = (a.decharge_lambdas
                            or getattr(fep_cfg, "decharge_lambdas",
                                       fep_cfg.lambdas))
                elif stage == "vdw":
                    lams = (a.vdw_lambdas
                            or getattr(fep_cfg, "vdw_lambdas",
                                       fep_cfg.lambdas))
                else:
                    lams = fep_cfg.lambdas
                reps = heal_leg(leg_dir, lams, rerun,
                                a.self_heal_attempts, logger=log)
                n_ok = sum(1 for r in reps.values() if r.completed)
                log.info("self-heal %s: %d/%d windows complete",
                         leg_name, n_ok, len(reps))
                for lam, r in reps.items():
                    if not r.completed:
                        last = r.history[-1] if r.history else {}
                        log.warning("  %s lambda_%s STILL FAILING (class=%s)",
                                    leg_name, lam,
                                    last.get("failure_class", "?"))

    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="amber-fep", description="ABFE/RBFE driver v2.4.21")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--submit",  action="store_true")
    p.add_argument("--skip-mdin-check", action="store_true",
                   help="Skip the v2.5.31 pre-submit mdin validator gate "
                        "(not recommended).")
    p.add_argument("--resume",  action="store_true",
                   help="Resubmit only unfinished windows (killed by "
                        "maintenance/preemption) and re-wire the "
                        "analyzer + cycle-closer. Skips topology "
                        "rebuild; reuses existing prmtops/windows.")
    p.add_argument("--analyze", action="store_true")
    p.add_argument("--mode", choices=["legacy","abfe","rbfe"], default="legacy")
    p.add_argument("--protein-pdb", type=Path)
    p.add_argument("--ligand-file","--ligand-sdf",dest="ligand_file",type=Path)
    p.add_argument("--ligand-index", type=int, default=0)
    p.add_argument("--ligand-name", type=str, default=None)
    p.add_argument("--ligand-charge", type=int, default=0)
    p.add_argument("--charge-method", default="bcc", choices=["bcc","gas","resp"])
    # v2.4.21: solvation box buffer (A). Threaded into SystemConfig for
    # both complex and solvent topology builds. Default 16.0 matches
    # SystemConfig and clears the pmemd.cuda small-box threshold.
    p.add_argument("--box-buffer", type=float, default=12.0,  # v2.5.71: 16->12 (see config.py rationale)
                   help="Solvation buffer in Angstrom (default 16.0).")
    p.add_argument("--lambdas", type=float, nargs="+",
                   default=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    # PATCH Bug 6 (v2.4.19): explicit per-stage schedules for ABFE.
    p.add_argument("--decharge-lambdas", type=float, nargs="+", default=None,
                   help="ABFE decharge-stage lambda schedule. If omitted, "
                        "FEPConfig.decharge_lambdas defaults are used.")
    p.add_argument("--vdw-lambdas", type=float, nargs="+", default=None,
                   help="ABFE vdW-stage lambda schedule. If omitted, "
                        "FEPConfig.vdw_lambdas defaults are used.")
    p.add_argument("--nstlim-eq",   type=int, default=250_000)
    p.add_argument("--nstlim-prod", type=int, default=10_000_000)  # v2.5.63: 10 ns @1fs true default
    p.add_argument("--restraint-leg", dest="restraint_leg",
                   action="store_true", default=False,
                   help="Run the explicit ~2 ns Boresch restraint MD leg and take its analytic dG from there (BAT.py-style validation). OFF by default.")
    p.add_argument("--no-restraint-leg", dest="restraint_leg",
                   action="store_false",
                   help="DEFAULT. Skip the ~0-FE restraint MD leg; fold the analytic Boresch standard-state term onto complex_vdw (Boresch-2003 / FEP+ / OpenFE style). The restraint POTENTIAL is still held ON during decharge+vdw.")
    p.add_argument("--restraint-nstlim-prod", type=int, default=2_000_000,
                   help="Option-A restraint leg production steps (default 2_000_000 = 2 ns @1fs; equilibration-only, not sampling-limited)")
    p.add_argument("--dt",          type=float, default=0.001)
    p.add_argument("--cutoff",      type=float, default=10.0)
    p.add_argument("--ntpr",        type=int, default=1000)
    p.add_argument("--ntwx",        type=int, default=5000)
    p.add_argument("--temperature", type=float, default=298.0)
    p.add_argument("--timask1", default=":LIG")
    p.add_argument("--timask2", default=":MOD")
    p.add_argument("--scmask1", default=":LIG")
    p.add_argument("--scmask2", default=":MOD")
    p.add_argument("--complex-prmtop",  type=Path)
    p.add_argument("--complex-inpcrd",  type=Path)
    p.add_argument("--solvent-prmtop",  type=Path)
    p.add_argument("--solvent-inpcrd",  type=Path)
    p.add_argument("--absolute-prmtop", type=Path)
    p.add_argument("--absolute-inpcrd", type=Path)
    p.add_argument("--boresch-json",    type=Path)
    p.add_argument("--auto-boresch",    action="store_true")
    p.add_argument("--fine-restraint", action="store_true",
                   help="Use front-loaded 20-window restraint schedule "
                        "(restraint_lambdas_fine) for clash-prone ligands.")
    p.add_argument("--auto-boresch-pdb", type=Path)
    p.add_argument("--ligand-resname",  default=None,
                   help="Ligand residue name used for ABFE masks and "
                        "tleap. If omitted, auto-detected from a .mol2 "
                        "(else defaults to LIG).")
    p.add_argument("--hremd", action="store_true")
    p.add_argument("--exchange-freq", type=int, default=1000)
    p.add_argument("--no-solvent", action="store_true")
    # ---- self-heal (final61): ensure every window produces a complete prod.out
    p.add_argument("--self-heal", action="store_true",
                   help="After analyze, ensure every lambda window of every "
                        "ABFE leg has a complete prod.out; diagnose failures, "
                        "adjust mdin parameters, and re-run.")
    p.add_argument("--self-heal-mode", choices=["bsub", "local", "dry-run"],
                   default="bsub",
                   help="How to re-run a window during self-heal.")
    p.add_argument("--no-auto-recover", action="store_true",
                   help="Disable the v2.5.15 completeness gate "
                        "that recovers incomplete windows "
                        "before analysis. Analysis then runs "
                        "on whatever windows exist.")
    p.add_argument("--self-heal-attempts", type=int, default=3,
                   help="Max remediation attempts per window (default 3).")
    # ---- RBFE (--mode rbfe) inputs (final58, Option B) ----
    p.add_argument("--ligands", type=Path, default=None,
                   help="RBFE: multi-record SDF/MOL2 of all ligands (same file "
                        "fed to rbfe_map).")
    p.add_argument("--rbfe-map-dir", type=Path, default=None,
                   help="RBFE: dir produced by amber_md.rbfe_map (edges.csv). "
                        "Defaults to <work-dir>/rbfe_map.")
    p.add_argument("--edges-csv", type=Path, default=None,
                   help="RBFE: explicit edges.csv path (overrides map-dir).")
    p.add_argument("--only-edges", nargs="*", default=None,
                   help="RBFE: restrict to these edge names (e.g. ligA~ligB).")
    p.add_argument("--project",   default=None,
                   help="LSF project (-P). Defaults to "
                        "HPCConfig.project (your-project).")
    p.add_argument("--queue",     default="gpu")
    p.add_argument("--queue-cpu", default=None)
    p.add_argument("--walltime",  default="24:00")
    p.add_argument("--n-gpu",     type=int, default=1)
    p.add_argument("--modules",   nargs="*", default=["amber/22"])
    p.add_argument("--venv",      default=None)
    a = p.parse_args(argv)
    if a.mode == "rbfe":
        return run_rbfe(a)
    return run_fep(a)


if __name__ == "__main__":
    sys.exit(main())
