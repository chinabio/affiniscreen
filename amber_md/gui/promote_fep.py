"""
promote_fep.py  --  Streamlit-free "promote MM-GBSA hit -> FEP" orchestration
(v2.5.0, Phase 4c).

Extracts the proven scaffold-and-submit logic from the legacy 8_Results.py
Promote-to-FEP section so the new Results pages can reuse one implementation.
Builds an amber_md.fep_driver command per ligand (relative TI or absolute
ABFE + Boresch), writes a submit_fep.sh, and runs the driver to scaffold.

Pure logic: no `import streamlit`.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import math
import shlex
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Lambda schedules
# ---------------------------------------------------------------------------
def linear_lambdas(n: int) -> list[float]:
    return [round(i / (n - 1), 3) for i in range(n)]


def cosine_lambdas(n: int) -> list[float]:
    """Denser at the endpoints (cosine spacing)."""
    return [round(0.5 * (1 - math.cos(math.pi * i / (n - 1))), 3)
            for i in range(n)]


def make_lambdas(n: int, spacing: str = "linear") -> list[float]:
    return cosine_lambdas(n) if spacing.startswith("denser") else linear_lambdas(n)


# ---------------------------------------------------------------------------
# Build-artifact discovery (prmtop/inpcrd/receptor) under a ligand dir
# ---------------------------------------------------------------------------
def find_fep_build(lig_dir: Path) -> dict:
    found: dict = {}
    build = Path(lig_dir) / "build"
    if build.is_dir():
        for kind in ("complex", "solvent"):
            for ext_prm, ext_crd in (("prmtop", "inpcrd"), ("parm7", "rst7")):
                prm = build / f"{kind}.{ext_prm}"
                crd = build / f"{kind}.{ext_crd}"
                if prm.exists() and crd.exists():
                    found[f"{kind}_prmtop"] = prm
                    found[f"{kind}_inpcrd"] = crd
                    break
        for cand in ("complex_solv.pdb", "complex.pdb", "receptor.pdb",
                     "protein.pdb"):
            p = build / cand
            if p.exists():
                found["receptor_pdb"] = p
                break
    return found


def has_fep_scaffold(wd: Path) -> bool:
    wd = Path(wd)
    return any((wd / f"fep_{wd.name}").glob("*")) if (wd / f"fep_{wd.name}").is_dir() else \
        any(wd.glob("fep_*/submit_fep.sh"))


# ---------------------------------------------------------------------------
# Build the fep_driver command for one ligand (no side effects).
# ---------------------------------------------------------------------------
def build_scaffold_cmd(lig_dir: Path, batch_dir: Path, opts: dict):
    """Return (cmd:list[str]|None, fep_wd:Path, reason:str). cmd is None when
    no usable build artifacts exist."""
    lig_dir, batch_dir = Path(lig_dir), Path(batch_dir)
    fep_wd = batch_dir / f"fep_{lig_dir.name}"
    build = find_fep_build(lig_dir)
    cmd = [sys.executable, "-m", "amber_md.fep_driver",
           "--work-dir", str(fep_wd),
           "--temperature", str(opts["temperature"]),
           "--nstlim-prod", str(int(opts["nstlim_prod"])),
           "--lambdas", *[f"{l:.3f}" for l in opts["lambdas"]]]
    if opts.get("hremd"):
        cmd += ["--hremd", "--exchange-freq", str(int(opts["exchange_freq"]))]

    have_legs = False
    if opts["mode"] == "relative":
        if "complex_prmtop" in build and "solvent_prmtop" in build:
            if opts.get("use_per_ligand_masks") and \
               opts.get("per_ligand_masks", {}).get(lig_dir.name):
                masks = opts["per_ligand_masks"][lig_dir.name]
            else:
                masks = {k: opts[k] for k in
                         ("timask1", "timask2", "scmask1", "scmask2")}
            cmd += ["--complex-prmtop", str(build["complex_prmtop"]),
                    "--complex-inpcrd", str(build["complex_inpcrd"]),
                    "--solvent-prmtop", str(build["solvent_prmtop"]),
                    "--solvent-inpcrd", str(build["solvent_inpcrd"]),
                    "--timask1", masks["timask1"],
                    "--timask2", masks["timask2"],
                    "--scmask1", masks["scmask1"],
                    "--scmask2", masks["scmask2"]]
            have_legs = True
    else:  # absolute
        if "complex_prmtop" in build:
            cmd += ["--mode", "abfe",
                    "--absolute-prmtop", str(build["complex_prmtop"]),
                    "--absolute-inpcrd", str(build["complex_inpcrd"]),
                    "--ligand-resname", opts["ligand_resname"],
                    "--auto-boresch"]
            if "receptor_pdb" in build:
                cmd += ["--auto-boresch-pdb", str(build["receptor_pdb"])]
            have_legs = True

    if not have_legs:
        return None, fep_wd, f"no usable prmtop/inpcrd in {lig_dir}/build/"

    cmd += ["--project", opts["project"], "--queue", opts["queue"],
            "--walltime", opts["walltime"]]
    if opts.get("modules"):
        cmd += ["--modules", *opts["modules"]]
    if opts.get("venv"):
        cmd += ["--venv", opts["venv"]]
    return cmd, fep_wd, ""


# ---------------------------------------------------------------------------
# Scaffold one ligand: write submit_fep.sh + run the driver. Returns a status
# dict identical in shape to the legacy _scaffold_fep.
# ---------------------------------------------------------------------------
def scaffold_fep(lig_dir: Path, batch_dir: Path, opts: dict,
                 timeout: int = 180, run_fn=None) -> dict:
    lig_dir = Path(lig_dir)
    cmd, fep_wd, reason = build_scaffold_cmd(lig_dir, batch_dir, opts)
    if cmd is None:
        return {"ligand": lig_dir.name, "status": "SKIP", "reason": reason,
                "cmd": ""}
    fep_wd.mkdir(parents=True, exist_ok=True)
    submit_sh = fep_wd / "submit_fep.sh"
    submit_sh.write_text(
        "#!/usr/bin/env bash\n"
        "# Auto-generated by Results 'Promote MM-GBSA -> FEP'.\n"
        "set -euo pipefail\n\n"
        + " \\\n  ".join(shlex.quote(c) for c in cmd) + " --submit\n")
    submit_sh.chmod(0o755)
    runner = run_fn or (lambda c: subprocess.run(
        c, capture_output=True, text=True, timeout=timeout))
    try:
        cp = runner(cmd)
    except subprocess.TimeoutExpired:
        return {"ligand": lig_dir.name, "status": "TIMEOUT",
                "reason": f"fep_driver scaffold > {timeout}s",
                "cmd": " ".join(shlex.quote(c) for c in cmd)}
    if cp.returncode != 0:
        out = (cp.stderr or cp.stdout) or ""
        return {"ligand": lig_dir.name, "status": "FAIL",
                "reason": out.splitlines()[-1][:240] if out else "(no output)",
                "cmd": " ".join(shlex.quote(c) for c in cmd)}
    return {"ligand": lig_dir.name, "status": "OK", "scaffold": str(fep_wd),
            "submit_sh": str(submit_sh),
            "cmd": " ".join(shlex.quote(c) for c in cmd)}
