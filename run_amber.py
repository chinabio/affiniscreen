#!/usr/bin/env python
"""AffiniScreen CLI - in-package analysis_kit / mmpbsa_report discovery."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import print_function
import sys, os, subprocess
from pathlib import Path

# Version resolved from the single source of truth (the VERSION file next to
# this script) so it can never drift from amber_md.__version__.
def _resolve_version():
    try:
        return (Path(__file__).resolve().parent / "VERSION").read_text().strip()
    except Exception:
        return "2.6.0"
__version__ = _resolve_version()

if sys.version_info < (3, 7):
    sys.stderr.write(
        "\n[amber_md] ERROR: needs Python >= 3.7 (found {}.{}.{}).\n"
        "Fix on the login node:\n"
        "    source activate_amber_md.sh\n"
        "Or first time:\n"
        "    bash bootstrap_fantasy.sh\n\n".format(
            sys.version_info[0], sys.version_info[1], sys.version_info[2]))
    sys.exit(2)

import argparse
HERE = Path(__file__).resolve().parent

def cmd_bootstrap():
    script = HERE/"bootstrap_fantasy.sh"
    if not script.exists(): sys.exit(f"ERROR: {script} not found")
    print(f"[run_amber] Running {script} ...")
    sys.exit(subprocess.call(["bash", str(script)]))

def cmd_activate():
    s = HERE/"activate_amber_md.sh"
    if not s.exists(): sys.exit(f"ERROR: {s} not found. Run --bootstrap first.")
    print(f"source {s}")

def main():
    p = argparse.ArgumentParser(
        description=f"AffiniScreen v{__version__} (Amber MD engine)",
        epilog="First time: python run_amber.py --bootstrap")
    p.add_argument("--version", action="version",
                   version=f"amber_md_workflow {__version__}")
    p.add_argument("--bootstrap", action="store_true",
                   help="Run bootstrap_fantasy.sh and exit.")
    p.add_argument("--activate", action="store_true",
                   help="Print the source command for activate_amber_md.sh and exit.")
    p.add_argument("--check-env", action="store_true",
                   help="Run env check and exit.")
    p.add_argument("--skip-env-check", action="store_true",
                   help="Skip pre-flight env check.")
    p.add_argument("--config", type=Path)
    p.add_argument("--pdb", type=Path)
    p.add_argument("--lig-resname", default=None,
                   help="Ligand residue name for MM-GBSA (forced via "
                        "antechamber -rn / PDB rewrite). If omitted, "
                        "auto-detected from a .mol2 ligand (else LIG).")
    p.add_argument("--protein-file", type=Path,
                   help="Protein-only file (.pdb or .mol2). Use with --ligand-file for DIRECT mode.")
    p.add_argument("--ligand-file", type=Path,
                   help="Ligand file (.sdf/.mol2/.pdb/.mol/.xyz). Format auto-detected.")
    p.add_argument("--ligand-charge", type=int, default=0,
                   help="Net ligand charge (default 0).")
    p.add_argument("--charge-method", default="bcc",
                   choices=["bcc","gas","resp"], help="Antechamber charge method.")
    p.add_argument("--ion-method", default="rand", choices=["rand","grid"],
                   help="Ion placement: 'rand' (fast) [default] or 'grid' (slow).")
    p.add_argument("--fast-ions", dest="ion_method", action="store_const",
                   const="rand", help="Shortcut for --ion-method=rand.")
    p.add_argument("--slow-ions", dest="ion_method", action="store_const",
                   const="grid", help="Shortcut for --ion-method=grid.")
    p.add_argument("--salt", type=float, default=0.15,
                   help="NaCl concentration in M (default 0.15).")
    p.add_argument("--no-protonation", dest="auto_protonation",
                   action="store_false", default=True,
                   help="Disable auto-detection of HIS/GLU/ASP protonation.")
    p.add_argument("--protonate", action="append", default=[],
                   metavar="CHAIN:RESNUM:NAME",
                   help="Manual override, e.g. --protonate D:198:GLH. Repeatable.")
    # ----- LSF / HPC overrides -----
    p.add_argument("--project", default=None,
                   help="LSF -P project (default: your-project).")
    p.add_argument("--queue", default=None,
                   help="LSF GPU queue (default: gpu).")
    p.add_argument("--walltime", default=None,
                   help="LSF walltime HH:MM (default: 24:00).")
    p.add_argument("--n-gpu", type=int, default=None,
                   help="GPUs per job (default: 1).")
    p.add_argument("--job-name", default=None,
                   help="LSF -J job name prefix (default: amberMD).")
    p.add_argument("--workdir", type=Path, default=Path("./run"))
    p.add_argument("--prod-ns", type=float, default=10.0)  # final51: unified 10 ns (both engines)
    p.add_argument("--equil-ns", type=float, default=1.0)
    p.add_argument("--no-submit",  action="store_true")
    p.add_argument("--no-monitor", dest="monitor",
                   action="store_false", default=True,
                   help="Submit the GPU job and exit; skip Stage 4/5. Default: monitor.")
    p.add_argument("--monitor", dest="monitor", action="store_true",
                   help="(Default) Wait for the GPU job and run Stage 4/5.")
    p.add_argument("--no-gbsa", action="store_true")
    p.add_argument("--decomp", action="store_true",
                   help="Enable MM/GBSA per-residue decomposition "
                        "(produces FINAL_DECOMP_MMPBSA.dat; adds ~30%% runtime).")
    p.add_argument("--decomp-residues", default="", metavar="MASK",
                   help="Amber residue mask, e.g. ':300-450', ':138,142,201'. "
                        "If empty with --decomp, MMPBSA defaults to all protein residues.")
    p.add_argument("--fep", action="store_true")
    p.add_argument("--fep-complex-prmtop", type=Path)
    p.add_argument("--fep-complex-inpcrd", type=Path)
    p.add_argument("--fep-solvent-prmtop", type=Path)
    p.add_argument("--fep-solvent-inpcrd", type=Path)
    p.add_argument("--fep-lambdas",
                   default="0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0")
    p.add_argument("--fep-timask1", default=":LIA")
    p.add_argument("--fep-timask2", default=":LIB")
    a = p.parse_args()

    proton_overrides = {}
    for spec in a.protonate:
        try:
            chain, resnum, name = spec.split(":")
            proton_overrides[(chain, int(resnum))] = name.upper()
        except ValueError:
            sys.exit(f"Bad --protonate value {spec!r}; expected CHAIN:RESNUM:NAME")
    if a.bootstrap: cmd_bootstrap()
    if a.activate: cmd_activate(); sys.exit(0)

    from amber_md.env import (check_python, check_python_packages,
                              check_amber_binaries, check_lsf, main as env_main)
    from amber_md.config import (WorkflowConfig, SystemConfig, MDConfig, MMGBSAConfig,
                                 FEPWorkflowConfig, FEPConfig, HPCConfig)
    from amber_md.pipeline import AmberPipeline

    if a.check_env: sys.exit(env_main([]))

    if not a.skip_env_check:
        print("[amber_md] Pre-flight environment check...")
        ok = check_python() and check_python_packages() and check_amber_binaries()
        check_lsf()
        if not ok:
            print("[amber_md] Pre-flight failed - try:")
            print("    source activate_amber_md.sh    # or --bootstrap first time")
            print("    python check_env.py            # full report")
            print("    --skip-env-check               # bypass")
            sys.exit(2)
        print("[amber_md] Pre-flight OK.\n")

    if a.config:
        cfg = WorkflowConfig.load(a.config)
        # final52: CLI MD-time flags override the loaded config when explicitly
        # passed, so --prod-ns / --equil-ns are never silently ignored on the
        # --config path (the GUI passes both). argparse defaults (10/1 ns) do
        # NOT override; only values that differ from the default are applied.
        if a.prod_ns is not None and abs(a.prod_ns - 10.0) > 1e-9:
            cfg.md.prod_nsteps = int(round(a.prod_ns * 1e6 / 2))
        if a.equil_ns is not None and abs(a.equil_ns - 1.0) > 1e-9:
            cfg.md.equil_nsteps = int(round(a.equil_ns * 1e6 / 2))
    else:
        has_combined = bool(a.pdb)
        has_direct   = bool(a.protein_file and a.ligand_file)
        if not (has_combined or has_direct):
            sys.exit("Need either --pdb, OR --protein-file + --ligand-file "
                     "(or --config / --bootstrap / --activate / --check-env)")
        fep = FEPWorkflowConfig(
            enabled=a.fep,
            params=FEPConfig(lambdas=tuple(float(x) for x in a.fep_lambdas.split(",")),
                             timask1=a.fep_timask1, timask2=a.fep_timask2,
                             scmask1=a.fep_timask1, scmask2=a.fep_timask2),
            complex_prmtop=a.fep_complex_prmtop, complex_inpcrd=a.fep_complex_inpcrd,
            solvent_prmtop=a.fep_solvent_prmtop, solvent_inpcrd=a.fep_solvent_inpcrd)
        _hpc_kwargs = {}
        if a.project   is not None: _hpc_kwargs["project"]   = a.project
        if a.queue     is not None: _hpc_kwargs["queue_gpu"] = a.queue
        if a.walltime  is not None: _hpc_kwargs["walltime"]  = a.walltime
        if a.n_gpu     is not None: _hpc_kwargs["n_gpu"]     = a.n_gpu
        if a.job_name  is not None: _hpc_kwargs["job_name"]  = a.job_name
        hpc = HPCConfig(**_hpc_kwargs)
        if hpc.venv_activate == "./activate_amber_md.sh":
            here_act = HERE/"activate_amber_md.sh"
            if here_act.exists():
                hpc.venv_activate = str(here_act)

        mmgbsa_cfg = MMGBSAConfig(
            enabled=not a.no_gbsa,
            decomposition=a.decomp,
            decomp_residues=a.decomp_residues or "")

        # Enforce the ligand residue name for MM-GBSA: explicit --lig-resname >
        # auto-detect from a .mol2 ligand > "LIG". WorkflowConfig/MMGBSA then
        # force it via antechamber -rn and a ligand-PDB rewrite, so prep is
        # self-consistent regardless of the residue name inside the input file.
        if not a.lig_resname:
            from amber_md.utils import detect_ligand_resname
            _lig = a.ligand_file if a.ligand_file else a.pdb
            if _lig:
                _rn, _src = detect_ligand_resname(_lig)
                a.lig_resname = _rn
                print("[run_amber] MM-GBSA ligand resname: %s (auto-detected from %s)"
                      % (_rn, _src))
            else:
                a.lig_resname = "LIG"
                print("[run_amber] MM-GBSA ligand resname: LIG (default)")
        else:
            print("[run_amber] MM-GBSA ligand resname: %s (user-specified)"
                  % a.lig_resname)
        cfg = WorkflowConfig(
            work_dir=a.workdir,
            complex_pdb=(a.pdb if a.pdb else Path("/dev/null")),
            ligand_resname=a.lig_resname,
            ligand_input=a.ligand_file, protein_input=a.protein_file,
            system=SystemConfig(ligand_charge=a.ligand_charge,
                                charge_method=a.charge_method,
                                ion_method=a.ion_method,
                                salt_conc_M=a.salt),
            md=MDConfig(prod_nsteps=int(a.prod_ns*1e6/2),
                        equil_nsteps=int(a.equil_ns*1e6/2)),
            hpc=hpc,
            mmgbsa=mmgbsa_cfg,
            fep=fep, submit=not a.no_submit, monitor=a.monitor,
            auto_protonation=a.auto_protonation,
            protonation_overrides=proton_overrides or None)
    AmberPipeline(cfg).run()

if __name__ == "__main__":
    main()
