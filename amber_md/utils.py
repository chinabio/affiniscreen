
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import shutil, subprocess, os, shlex
from pathlib import Path
from .logger import get_logger
log = get_logger()
class CommandError(RuntimeError): pass
def which_or_die(t):
    p = shutil.which(t)
    if not p: raise CommandError(f"Required executable '{t}' not found. Did you `module load amber`?")
    return p
def ensure_dir(p): p = Path(p); p.mkdir(parents=True, exist_ok=True); return p
def run(cmd, cwd=None, stdin=None, env=None, check=True, capture=True):
    cmd_list = shlex.split(cmd) if isinstance(cmd, str) else list(cmd)
    log.info("RUN: %s%s", " ".join(cmd_list), f"  (cwd={cwd})" if cwd else "")
    try:
        cp = subprocess.run(cmd_list, cwd=str(cwd) if cwd else None, input=stdin,
                            text=True, capture_output=capture,
                            env={**os.environ, **(env or {})}, check=False)
    except FileNotFoundError as e:
        raise CommandError(str(e)) from e
    if check and cp.returncode != 0:
        raise CommandError(f"Command failed (exit {cp.returncode}): {' '.join(cmd_list)}\nSTDERR:\n{cp.stderr or ''}")
    return cp


def detect_ligand_resname(ligand_file, index: int = 0, default: str = "LIG"):
    """Best-effort ligand residue name read from the input file.

    * .mol2 -> residue-name column (8th field) of the @<TRIPOS>ATOM block,
      e.g. 'UNK'. Trailing digits stripped; capped at 4 chars (PDB/tleap).
    * .sdf/.mol -> SDF carries no residue name; antechamber/tleap assigns one,
      so the default is returned.
    Returns (resname, source) with source in {'mol2','default'}. Never raises.
    Shared single implementation used by the CLI driver and the GUI.
    """
    from pathlib import Path as _Path
    try:
        p = _Path(ligand_file); ext = p.suffix.lower()
        if ext == ".mol2":
            in_atoms = False
            for line in p.read_text(errors="ignore").splitlines():
                s = line.strip()
                if s.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True; continue
                if in_atoms:
                    if s.startswith("@<TRIPOS>"):
                        break
                    cols = s.split()
                    if len(cols) >= 8:
                        rn = cols[7].strip()
                        rn = "".join(ch for ch in rn if not ch.isdigit()) or rn
                        if rn:
                            return rn[:4], "mol2"
        return default, "default"
    except Exception:
        return default, "default"
