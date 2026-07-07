
# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

# -*- coding: utf-8 -*-
"""Environment checker / setup helper.

Python 2.7 / 3.5+ compatible so it can DIAGNOSE old Python before any
3.7+-only imports (dataclasses) are triggered.

Run:
    python -m amber_md.env
    python check_env.py
"""
from __future__ import print_function
import os, sys, shutil, subprocess

PYTHON_MIN = (3, 7)

REQ_PY_PACKAGES = [
    ("numpy", "numpy", "RMSD/RMSF + FEP integration"),
]
OPT_PY_PACKAGES = [
    ("matplotlib", "matplotlib", "PNG plots"),
    ("parmed",     "parmed",     "Topology split (provided by amber/22.8)"),
    ("MDAnalysis", "MDAnalysis", "Viz fallback (b)"),
    ("mdtraj",     "mdtraj",     "NGLview HTML (b)"),
    ("nglview",    "nglview",    "Interactive HTML (b)"),
    ("pymol",      None,         "PNG snapshots — module load pymol/3.0.4"),
    ("streamlit",  "streamlit",  "GUI"),
    ("pymbar",     "pymbar",     "BAR estimator for FEP (TI works without)"),
]
REQ_AMBER_BINS = [
    ("tleap",        "Build solvated system"),
    ("antechamber",  "Ligand parametrization"),
    ("parmchk2",     "Ligand frcmod"),
    ("pdb4amber",    "PDB cleaning"),
    ("pmemd.cuda",   "GPU MD engine"),
    ("cpptraj",      "RMSD/RMSF analysis"),
    ("MMPBSA.py",    "MM/GBSA"),
]
OPT_AMBER_BINS = [
    ("ante-MMPBSA.py", "Topology split fallback when ParmEd missing"),
]
REQ_LSF_BINS = [("bsub","Job submission"), ("bjobs","Job monitoring")]

GREEN  = "\033[32m" if sys.stdout.isatty() else ""
RED    = "\033[31m" if sys.stdout.isatty() else ""
YELLOW = "\033[33m" if sys.stdout.isatty() else ""
RESET  = "\033[0m"  if sys.stdout.isatty() else ""

def ok(m):    print("  " + GREEN  + "[OK]  " + RESET + m)
def bad(m):   print("  " + RED    + "[FAIL]" + RESET + " " + m)
def warn(m):  print("  " + YELLOW + "[WARN]" + RESET + " " + m)
def info(m):  print("        " + m)

def have_module(name):
    try:
        import importlib
        if hasattr(importlib, "util") and hasattr(importlib.util, "find_spec"):
            return importlib.util.find_spec(name) is not None
        import imp
        try: imp.find_module(name); return True
        except ImportError: return False
    except Exception: return False

def have_bin(name):
    if hasattr(shutil, "which"):
        return shutil.which(name)
    for d in os.environ.get("PATH","").split(os.pathsep):
        cand = os.path.join(d, name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None

def check_python():
    print("\n== Python interpreter ==")
    v = sys.version_info; full = ".".join(str(x) for x in v[:3])
    if (v.major, v.minor) >= PYTHON_MIN:
        ok("Python {} at {}".format(full, sys.executable)); return True
    bad("Python {} too old — need >= {}.{}".format(full, PYTHON_MIN[0], PYTHON_MIN[1]))
    info("This pipeline uses 'dataclasses' (stdlib in Python 3.7+).")
    info("Fix on a cluster:")
    info("  module avail python")
    info("  module load python/3.11    # or whichever 3.7+ module exists")
    info("Or use Amber's bundled Python:")
    info("  module load amber && which python && python --version")
    info("Or with conda:")
    info("  conda create -n amber-md -c conda-forge python=3.11")
    info("  conda activate amber-md")
    return False

def check_python_packages():
    print("\n== Python packages ==")
    all_ok = True
    for mod, pip, why in REQ_PY_PACKAGES:
        if have_module(mod): ok("{:12s} ({})".format(mod, why))
        else:
            bad("{:12s} MISSING — {}".format(mod, why))
            info("Fix:  pip install {}  (or use bootstrap_fantasy.sh)".format(pip)); all_ok = False
    for mod, pip, why in OPT_PY_PACKAGES:
        if have_module(mod): ok("{:12s} ({})".format(mod, why))
        else:
            install = "pip install " + pip if pip else "module load pymol (or conda)"
            warn("{:12s} optional — {}.  {}".format(mod, why, install))
    return all_ok

def check_amber_binaries():
    print("\n== Amber binaries ==")
    all_ok = True
    for b, why in REQ_AMBER_BINS:
        p = have_bin(b)
        if p: ok("{:14s} {}  ({})".format(b, p, why))
        else:
            bad("{:14s} not on PATH — {}".format(b, why)); all_ok = False
    if not all_ok:
        info("Fix on the login node:")
        info("  source activate_amber_md.sh")
        info("Or:  module load amber/22.8")
    for b, why in OPT_AMBER_BINS:
        p = have_bin(b)
        if p: ok("{:14s} {}  ({})".format(b, p, why))
        else: warn("{:14s} optional — {}".format(b, why))
    return all_ok

def check_lsf():
    print("\n== LSF (submission) ==")
    any_missing = False
    for b, why in REQ_LSF_BINS:
        p = have_bin(b)
        if p: ok("{:14s} {}  ({})".format(b, p, why))
        else:
            warn("{:14s} missing — {}. Use --no-submit to bypass.".format(b, why))
            any_missing = True
    return not any_missing

def check_gpu():
    print("\n== GPU (CUDA) ==")
    nvsmi = have_bin("nvidia-smi")
    if not nvsmi:
        warn("nvidia-smi not found — login node? GPU queue submission may still work.")
        return True
    try:
        out = subprocess.check_output(
            [nvsmi, "--query-gpu=name,memory.total,driver_version",
             "--format=csv,noheader"], stderr=subprocess.STDOUT)
        out = out.decode("utf-8", errors="replace").strip()
        if out:
            for line in out.splitlines(): ok("GPU: " + line.strip())
        else:
            warn("nvidia-smi returned no GPUs (login node?).")
    except Exception as e:
        warn("nvidia-smi failed: {}".format(e))
    return True

def check_env_vars():
    print("\n== Relevant environment variables ==")
    for var in ("AMBERHOME","CUDA_HOME","LD_LIBRARY_PATH","CUDA_VISIBLE_DEVICES"):
        val = os.environ.get(var)
        if val: ok("{:20s} = {}".format(var, val))
        else:   warn("{:20s} not set".format(var))
    return True

def write_setup_template(path):
    content = """#!/usr/bin/env bash
# Auto-generated by amber_md.env — edit module names for your cluster.
module purge
module load gcc/11.5
module load cuda
module load amber
module load python/3.11
if [ -f "$HOME/venvs/amber-md/bin/activate" ]; then
    source "$HOME/venvs/amber-md/bin/activate"
fi
python -m amber_md.env
"""
    with open(path,"w") as f: f.write(content)
    try: os.chmod(path, 0o755)
    except Exception: pass

def check_conda_env():
    """v2.2: check we are in the right conda env."""
    print("\n== Conda environment ==")
    cenv = os.environ.get("CONDA_DEFAULT_ENV")
    if cenv == "amber-md":
        ok("Conda env 'amber-md' is active")
    elif cenv:
        warn("Conda env '{}' active — expected 'amber-md'.".format(cenv))
        info("Fix: source activate_amber_md.sh")
    else:
        warn("No conda env active.")
        info("Fix: source activate_amber_md.sh")
    return True

def main(argv=None):
    argv = argv or sys.argv[1:]
    write_template = "--write-template" in argv
    print("=" * 64)
    print(" AffiniScreen — environment check")
    print("=" * 64)
    results = [("Python", check_python())]
    if results[-1][1]:
        results.append(("Python packages", check_python_packages()))
    results.append(("Amber binaries", check_amber_binaries()))
    results.append(("LSF",            check_lsf()))
    check_gpu(); check_env_vars(); check_conda_env()
    print("\n" + "=" * 64)
    failed = [n for n,ok_ in results if not ok_]
    if failed:
        print(RED + "FAILED:" + RESET + " " + ", ".join(failed))
        print("Re-run after fixing the items above.")
        print("the login node quick-fix:  source activate_amber_md.sh")
        rc = 1
    else:
        print(GREEN + "All required checks passed." + RESET); rc = 0
    if write_template:
        write_setup_template("setup_env.sh")
        print("\nWrote setup_env.sh — edit and `source` it before running.")
    print("=" * 64)
    return rc

if __name__ == "__main__":
    sys.exit(main())
