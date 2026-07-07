"""OpenFE+OpenMM helpers for the RBFE / ABFE Streamlit pages.
import logging as _logging
log = _logging.getLogger(__name__)

Lives at amber_md/gui/openfe_common.py.

Verified against the openfe 1.11.x CLI:
  * `plan-rbfe-network -s settings.yaml` configures ONLY mapper / network /
    partial_charge -- NOT simulation length, platform, or sampler. Those live in
    the Python protocol layer, so sim-time / platform overrides are applied by a
    generated planner script (amber_md/openfe_plan.py) instead of the YAML.
  * `--n-protocol-repeats=1` lets each repeat be submitted as its OWN quickrun
    job on the same edge JSON. We exploit this to parallelise repeats across the
    8-GPU `gpu` queue and keep per-job walltime short.
  * `gather --report dg -o out.tsv` writes a tab-separated absolute-deltaG table.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import streamlit as st

# v2.5.0: client-side GPU throttling removed -- the LSF queue handles job
# control. MAX_GPU retained only for backward-compat imports; 0 = no cap.
MAX_GPU = 0  # 0 = no client-side cap


@dataclass
class OpenFESettings:
    # protocol (applied via the Python planner, NOT the CLI YAML)
    sim_time_ns: float
    equil_time_ns: float
    n_replicates: int          # how many quickrun jobs we submit per edge
    sampler: str               # "repex" | "sams" | "independent"
    platform: str              # "CUDA" | "OpenCL"
    forcefield: str            # small-molecule FF, e.g. "openff-2.2.0"
    # mapper / network / charges (applied via the CLI -s YAML)
    mapper: str                # "KartografAtomMapper" | "LomapAtomMapper"
    network: str               # "generate_minimal_spanning_network" | ...
    charge_method: str         # "am1bcc" | "nagl" | "am1bccelf10" | "espaloma"
    # LSF
    project: str
    queue: str
    walltime: str
    conda_env: str
    extra_modules: list[str] = field(default_factory=list)
    max_concurrent: int = 0  # 0 = unlimited (LSF queues the rest)


def render_openfe_settings(key_prefix: str) -> OpenFESettings:
    st.markdown("#### Sampling (applied via the Python protocol layer)")
    c1, c2, c3 = st.columns(3)
    sim_time = c1.number_input(
        "Production / lambda (ns)", 0.1, 20.0, 5.0, 0.1, key=f"{key_prefix}_simns",
        help="OpenFE default 5 ns. Use <=0.5 ns for a smoke test.")
    equil_time = c2.number_input(
        "Equilibration (ns)", 0.0, 5.0, 1.0, 0.1, key=f"{key_prefix}_eqns")
    n_rep = c3.number_input(
        "Replicates (submitted as parallel jobs)", 1, 5, 3, 1,
        key=f"{key_prefix}_reps",
        help="Plan uses --n-protocol-repeats=1; we submit this many quickrun "
             "jobs per edge so repeats run in parallel across GPUs.")

    c1, c2, c3 = st.columns(3)
    sampler = c1.selectbox("Sampler", ["repex", "sams", "independent"], 0,
                           key=f"{key_prefix}_sampler",
                           help="repex = HREX (best convergence).")
    platform = c2.selectbox("Platform", ["CUDA", "OpenCL"], 0,
                            key=f"{key_prefix}_plat")
    forcefield = c3.selectbox(
        "Small-molecule FF",
        ["openff-2.2.0", "openff-2.1.0", "openff-2.0.0", "gaff-2.11"], 0,
        key=f"{key_prefix}_ff")

    st.markdown("#### Network planning (applied via the CLI `-s` YAML)")
    c1, c2, c3 = st.columns(3)
    mapper = c1.selectbox("Atom mapper",
                          ["KartografAtomMapper", "LomapAtomMapper"], 0,
                          key=f"{key_prefix}_mapper")
    network = c2.selectbox(
        "Network algorithm",
        ["generate_minimal_spanning_network",
         "generate_minimal_redundant_network",
         "generate_lomap_network",
         "generate_radial_network"], 0, key=f"{key_prefix}_net",
        help="minimal_redundant adds cycles -> better cycle-closure QC.")
    # nagl is available in your env (openff-nagl installed) -> fast ML charges
    charge_method = c3.selectbox(
        "Partial charges", ["am1bcc", "nagl", "am1bccelf10", "espaloma"], 0,
        key=f"{key_prefix}_chg",
        help="am1bcc (ambertools, default). nagl = fast ML charges (installed).")

    st.divider()
    st.markdown("#### HPC (LSF)")
    c1, c2 = st.columns(2)
    project = c1.text_input("LSF project (-P)", "your-project",
                            key=f"{key_prefix}_proj")
    queue = c2.text_input("GPU queue (-q)", "gpu", key=f"{key_prefix}_q")
    walltime = c1.text_input("Walltime (HH:MM)", "24:00", key=f"{key_prefix}_wt")
    # v2.5.0: no client-side throttle. The LSF queue dispatches as GPUs free
    # up; we submit everything at once. max_concurrent kept (=0) so code that
    # reads it still works.
    max_conc = 0
    conda_env = st.text_input("conda env", "openfe_env", key=f"{key_prefix}_env")
    extra_mods = st.text_input("extra `module load`", "cuda",
                               key=f"{key_prefix}_mods").split()

    if sim_time <= 0.6 and n_rep == 1:
        st.info("Smoke-test settings -- validates the pipeline, NOT for "
                "production numbers.")

    return OpenFESettings(
        sim_time_ns=float(sim_time), equil_time_ns=float(equil_time),
        n_replicates=int(n_rep), sampler=sampler, platform=platform,
        forcefield=forcefield, mapper=mapper, network=network,
        charge_method=charge_method, project=project, queue=queue,
        walltime=walltime, conda_env=conda_env, extra_modules=extra_mods,
        max_concurrent=int(max_conc))


# ---------------------------------------------------------------------------
# CLI `-s` YAML: ONLY mapper / network / partial_charge (verified scope)
# ---------------------------------------------------------------------------
def write_network_yaml(s: OpenFESettings, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(f"""
# Auto-generated. Scope per `openfe plan-rbfe-network -s`: mapper/network/charge ONLY.
mapper:
  method: {s.mapper}
network:
  method: {s.network}
partial_charge:
  method: {s.charge_method}
""")
    return out_path


def plan_rbfe_cmd(protein_pdb: Path, ligand_sdf: Path, out_dir: Path,
                  settings_yaml: Path | None, n_cores: int = 1,
                  openfe_bin: str = "openfe") -> list[str]:
    """Plan with --n-protocol-repeats=1 so repeats can be parallelised.

    openfe_bin lets callers point at a specific OpenFE executable
    (e.g. one inside a dedicated conda env) instead of relying on PATH.
    """
    cmd = [
        openfe_bin, "plan-rbfe-network",
        "-M", str(ligand_sdf),
        "-p", str(protein_pdb),
        "-o", str(out_dir),
        "--n-protocol-repeats", "1",   # KEY: one repeat per quickrun call
    ]
    if n_cores > 1:
        cmd += ["-n", str(n_cores)]
    if settings_yaml is not None:
        cmd += ["-s", str(settings_yaml)]
    return cmd


def quickrun_cmd(edge_json: Path, result_json: Path, work_dir: Path,
                 resume: bool = False) -> list[str]:
    cmd = ["openfe", "quickrun", str(edge_json),
           "-o", str(result_json), "-d", str(work_dir)]
    if resume:
        cmd.append("--resume")
    return cmd


def gather_cmd(results_dir: Path, out_tsv: Path,
               report: str = "dg") -> list[str]:
    """gather --report dg -o out.tsv --allow-partial (verified flags)."""
    return ["openfe", "gather", str(results_dir),
            "--report", report, "-o", str(out_tsv), "--allow-partial"]


# ---------------------------------------------------------------------------
# Per-(edge, repeat) LSF submission script
# ---------------------------------------------------------------------------
def make_edge_bsub_script(edge_json: Path, result_json: Path, work_dir: Path,
                          s: OpenFESettings, job_name: str,
                          log_dir: Path) -> str:
    log_dir.mkdir(parents=True, exist_ok=True)
    module_lines = "\n".join(f"module load {m}" for m in s.extra_modules if m)
    # FIX: do NOT --resume by default. quickrun refuses to run if result_json
    # already exists ("Path ... is a file"); a fresh dir + no-resume avoids
    # reusing a poisoned/partial previous attempt. The script wipes the work
    # dir below so every submission starts clean.
    qr = " ".join(str(c) for c in quickrun_cmd(edge_json, result_json, work_dir,
                                               resume=False))
    # Only emit an LSF project line when one is set; an empty -P is malformed
    # and gets rejected by the scheduler. (Every the login node bsub needs a real -P.)
    _proj = str(getattr(s, "project", "") or "").strip()
    _p_line = f"#BSUB -P {_proj}\n" if _proj else ""
    if not _p_line:
        log.warning("OpenFE submit: no LSF project set; omitting -P. The job may "
                    "be rejected -- set project (e.g. your-project).")
    return f"""#!/bin/bash
#BSUB -q {s.queue}
{_p_line}#BSUB -J {job_name}
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=shared:j_exclusive=yes"
#BSUB -W {s.walltime}
#BSUB -o {log_dir}/%J.out
#BSUB -e {log_dir}/%J.err

set -euo pipefail

module purge
{module_lines}

source ~/miniforge3/etc/profile.d/conda.sh
conda activate {s.conda_env}

echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${{CUDA_VISIBLE_DEVICES:-<unset>}}"
nvidia-smi --query-gpu=name,memory.total,compute_mode --format=csv || true

# Fail fast if OpenMM cannot see a CUDA device (instead of burning wall time)
python -c "import openmm; print('OpenMM', openmm.__version__); print('CUDA platform:', openmm.Platform.getPlatformByName('CUDA').getName())"

# FIX: fresh working dir so quickrun never collides with an existing results.json
rm -rf "{work_dir}"
mkdir -p "{work_dir}"

{qr}
echo "EDGE-REPEAT DONE: {job_name}"
"""


def count_my_jobs(job_name_prefix: str) -> int:
    user = os.environ.get("USER", "")
    try:
        cp = subprocess.run(
            ["bjobs", "-u", user, "-J", f"{job_name_prefix}*",
             "-o", "JOBID", "-noheader"],
            capture_output=True, text=True, timeout=10)
        return 0 if cp.returncode != 0 else \
            len([l for l in cp.stdout.splitlines() if l.strip()])
    except Exception:
        return 0


def submit_bsub(script_path: Path) -> tuple[bool, str]:
    try:
        with open(script_path) as f:
            cp = subprocess.run(["bsub"], stdin=f, capture_output=True,
                                text=True, timeout=30)
        return cp.returncode == 0, (cp.stdout or cp.stderr).strip()
    except Exception as e:
        return False, str(e)


def list_transformations(network_dir: Path) -> list[Path]:
    tdir = network_dir / "transformations"
    if tdir.is_dir():
        return sorted(tdir.glob("*.json"))
    return sorted(network_dir.glob("*.json"))


def edge_repeat_status(result_json: Path, work_dir: Path) -> str:
    if result_json.exists():
        try:
            data = json.loads(result_json.read_text())
            if data.get("estimate") is not None or "protocol_result" in data:
                return "DONE"
        except Exception:
            return "analyzing"
        return "analyzing"
    return "running" if work_dir.exists() else "not started"


# ---------------------------------------------------------------------------
# Preflight: AM1BCC charge generation needs a self-consistent AmberTools
# (antechamber) reachable from the SAME env that runs `openfe`. A common
# failure is the conda env shelling out to a broken system Amber, e.g.
# /share/apps/.../wrapped_progs/antechamber -> "AmberTools version None".
# This check catches that BEFORE we spawn the planner.
# ---------------------------------------------------------------------------
def preflight_openfe_charges(openfe_bin: str = "openfe",
                             charge_method: str = "am1bcc") -> dict:
    """Validate the env can generate the requested partial charges.

    Returns {"ok": bool, "errors": [str], "warnings": [str], "info": {str:str}}.

    IMPORTANT: this evaluates the *sanitized* environment that the GUI will
    actually spawn `openfe` with (see sanitized_openfe_env) -- NOT the GUI
    process's own (possibly polluted) os.environ. That way a system Amber that
    is merely inherited by the GUI does not block a launch that will be cleaned
    up at spawn time. If the GUI's raw env is dirty we note it informationally.

    Only AM1BCC-family methods require AmberTools/antechamber. RDKit methods
    (gasteiger/mmff94) and NAGL skip the antechamber checks.
    """
    import os, shutil, subprocess, sys, json
    from pathlib import Path

    errors: list[str] = []
    warnings: list[str] = []
    info: dict[str, str] = {}

    needs_amber = str(charge_method).lower() in {
        "am1bcc", "am1bccelf10", "am1-mulliken"}
    info["charge_method"] = str(charge_method)
    if not needs_amber:
        info["note"] = (f"charge_method '{charge_method}' does not require "
                        "AmberTools; skipping antechamber checks.")
        return {"ok": True, "errors": errors, "warnings": warnings, "info": info}

    ofe = shutil.which(openfe_bin) or openfe_bin
    info["openfe_bin"] = ofe
    ofe_path = Path(ofe)
    env_bin = ofe_path.parent                       # .../<env>/bin
    env_root = env_bin.parent
    env_python = ofe_path.with_name("python")
    if not env_python.exists():
        env_python = Path(sys.executable)
        warnings.append(
            f"Could not find a python next to '{ofe}'; using {env_python} for "
            "the toolkit probe (may not reflect the planner's env).")
    info["probe_python"] = str(env_python)

    # ---- Build the SANITIZED env: this is what the GUI actually launches with.
    sanitized = None
    try:
        from amber_md.gui.common import sanitized_openfe_env as _san
        sanitized = _san(openfe_bin)
    except Exception as e:  # pragma: no cover - fall back to raw env
        warnings.append(f"Could not build sanitized env ({e}); evaluating raw "
                        "os.environ instead.")
        sanitized = dict(os.environ)

    eff_path = sanitized.get("PATH", "")
    eff_pp = sanitized.get("PYTHONPATH", "")
    eff_amberhome = sanitized.get("AMBERHOME", "")

    # Informational: was the GUI's OWN env dirty? (the spawn fixes this)
    def _is_system_amber_dir(p: str) -> bool:
        from pathlib import Path as _P
        if not p:
            return False
        try:
            rp = _P(p).resolve()
        except Exception:
            rp = _P(p)
        s = str(rp).lower()
        if str(env_root.resolve()).lower() in s:
            return False
        parts = set(rp.parts)
        amber_root = any(
            seg.startswith("amber") and "_" not in seg and "-" not in seg
            and (seg[5:].isdigit() or seg in ("amber", "ambertools"))
            for seg in parts)
        if not amber_root:
            return False
        return ("wrapped_progs" in parts) or ("bin" in parts) or \
               ("site-packages" in parts and "lib" in parts)

    def _first_sys_amber(varval):
        for p in (varval or "").split(os.pathsep):
            if _is_system_amber_dir(p):
                return p
        return None

    raw_dirty = (_first_sys_amber(os.environ.get("PATH", "")) or
                 _first_sys_amber(os.environ.get("PYTHONPATH", "")))
    if raw_dirty:
        warnings.append(
            "The GUI process inherited a system Amber (e.g. an auto-loaded "
            f"`amber` module): '{raw_dirty}'. The GUI sanitizes this before "
            "spawning OpenFE, so the launch is OK -- but for a permanent fix "
            "`module unload amber` (and `module save`) so future logins/GUI "
            "starts are clean.")

    # ---- (a) antechamber as seen in the SANITIZED PATH.
    ac = shutil.which("antechamber", path=eff_path)
    info["antechamber (sanitized PATH)"] = ac or "<not found>"
    if ac is None:
        errors.append(
            "antechamber not found even after sanitizing PATH. Install "
            "AmberTools into the OpenFE env:  conda install -n <openfe_env> "
            "-c conda-forge ambertools")
    elif Path(ac).resolve().parent != env_bin.resolve():
        # Sanitized PATH still resolves antechamber outside the env -> real problem.
        errors.append(
            f"After sanitizing, antechamber still resolves to '{ac}', outside "
            f"the OpenFE env bin '{env_bin}'. Install ambertools into the env.")

    # ---- (b) AMBERHOME in the sanitized env.
    info["AMBERHOME (sanitized)"] = eff_amberhome or "<unset>"
    if eff_amberhome:
        try:
            inside = Path(eff_amberhome).resolve() == env_root.resolve()
        except Exception:
            inside = False
        if not inside:
            warnings.append(
                f"Sanitized AMBERHOME='{eff_amberhome}' is not the OpenFE env "
                f"('{env_root}').")

    # ---- (b2) system Amber remaining in the SANITIZED PYTHONPATH (must be none).
    info["PYTHONPATH (sanitized)"] = eff_pp or "<unset>"
    leftover_pp = _first_sys_amber(eff_pp)
    if leftover_pp:
        errors.append(
            f"A system Amber is STILL on PYTHONPATH after sanitizing: "
            f"'{leftover_pp}'. OpenFF would import parmed/antechamber from there "
            "and crash with 'AmberTools version None'.")

    # ---- (c) Authoritative probe under the SANITIZED env.
    probe = (
        "import json\n"
        "try:\n"
        "    from openff.toolkit.utils import AmberToolsToolkitWrapper as A\n"
        "    w=A(); print(json.dumps({'available':bool(w.is_available()),"
        "'version':getattr(w,'toolkit_version',None)}))\n"
        "except Exception as e:\n"
        "    print(json.dumps({'available':False,'error':f'{type(e).__name__}: {e}'}))\n"
    )
    try:
        cp = subprocess.run([str(env_python), "-c", probe],
                            capture_output=True, text=True, timeout=120,
                            env=sanitized)
        payload = {}
        for ln in (cp.stdout or "").splitlines():
            ln = ln.strip()
            if ln.startswith("{"):
                try:
                    payload = json.loads(ln)
                except Exception:
                    pass
        if payload.get("available"):
            info["ambertools_toolkit (sanitized)"] = (
                f"available (version {payload.get('version') or '?'})")
        else:
            err = payload.get("error", "is_available() returned False")
            errors.append(
                "Under the sanitized env, OpenFF AmberToolsToolkitWrapper is NOT "
                f"available -- AM1BCC would fail. Detail: {err}")
            if cp.stderr.strip():
                info["probe_stderr"] = cp.stderr.strip()[-800:]
    except Exception as e:  # noqa: BLE001
        warnings.append(
            f"Could not run the AmberTools toolkit probe ({type(e).__name__}: "
            f"{e}); proceeding without it.")

    return {"ok": not errors, "errors": errors,
            "warnings": warnings, "info": info}


def preflight_openmm_mmgbsa(probe_python: str = "python",
                            platform: str = "CUDA") -> dict:
    """Validate the env can run OpenMM-MD + AmberTools MMPBSA.py MM-GBSA.

    Mirrors preflight_openfe_charges. Checks, against the SAME python that will
    run the pipeline (the OpenFE env python, where OpenMM + AmberTools live):
      (a) AmberTools binaries on PATH (tleap/antechamber/parmchk2/pdb4amber/
          MMPBSA.py/ante-MMPBSA.py),
      (b) OpenMM importable,
      (c) the requested OpenMM Platform (e.g. CUDA) is available.

    Returns {"ok": bool, "errors": [str], "warnings": [str], "info": {str:str}}.
    """
    import os, shutil, subprocess, sys, json
    from pathlib import Path

    errors: list[str] = []
    warnings: list[str] = []
    info: dict[str, str] = {}
    info["requested_platform"] = platform

    py = shutil.which(probe_python) or probe_python
    if not Path(py).exists():
        py = sys.executable
        warnings.append(f"probe python '{probe_python}' not found; using {py}.")
    info["probe_python"] = py

    # Build the sanitized env so the probe reflects the actual launch env.
    sanitized = None
    try:
        from amber_md.gui.common import sanitized_openfe_env as _san
        sanitized = _san(py)
    except Exception:
        sanitized = dict(os.environ)
    eff_path = sanitized.get("PATH", os.environ.get("PATH", ""))

    # (a) AmberTools binaries.
    req = ["tleap", "antechamber", "parmchk2", "pdb4amber",
           "MMPBSA.py", "ante-MMPBSA.py"]
    missing = [t for t in req if shutil.which(t, path=eff_path) is None]
    info["ambertools"] = "all present" if not missing else f"MISSING: {missing}"
    if missing:
        errors.append(
            "AmberTools binaries missing from PATH: " + ", ".join(missing) +
            ". Install AmberTools into the env: conda install -n <env> "
            "-c conda-forge ambertools")

    # (b)+(c) OpenMM import + platform availability via the probe python.
    probe = (
        "import json\n"
        "out={}\n"
        "try:\n"
        "    import openmm\n"
        "    out['openmm_version']=openmm.__version__\n"
        "    plats=[openmm.Platform.getPlatform(i).getName() "
        "for i in range(openmm.Platform.getNumPlatforms())]\n"
        "    out['platforms']=plats\n"
        "    rep=None\n"
        "    try:\n"
        "        from openmm.app import NetCDFReporter as _R\n"
        "        rep='openmm.app.NetCDFReporter'\n"
        "    except Exception:\n"
        "        try:\n"
        "            from mdtraj.reporters import NetCDFReporter as _R\n"
        "            rep='mdtraj.reporters.NetCDFReporter'\n"
        "        except Exception:\n"
        "            rep=None\n"
        "    out['netcdf_reporter']=rep\n"
        "    out['ok']=True\n"
        "except Exception as e:\n"
        "    out={'ok':False,'error':f'{type(e).__name__}: {e}'}\n"
        "print(json.dumps(out))\n"
    )
    try:
        cp = subprocess.run([py, "-c", probe], capture_output=True, text=True,
                            timeout=120, env=sanitized)
        payload = {}
        for ln in (cp.stdout or "").splitlines():
            ln = ln.strip()
            if ln.startswith("{"):
                try:
                    payload = json.loads(ln)
                except Exception:
                    pass
        if payload.get("ok"):
            info["openmm_version"] = payload.get("openmm_version", "?")
            plats = payload.get("platforms", [])
            info["openmm_platforms"] = ", ".join(plats) or "<none>"
            if platform not in plats:
                errors.append(
                    f"Requested OpenMM platform '{platform}' is NOT available. "
                    f"Available: {plats}. Install the CUDA-enabled OpenMM build "
                    "or pick an available platform.")
            rep = payload.get("netcdf_reporter")
            info["netcdf_reporter"] = rep or "<none>"
            if not rep:
                errors.append(
                    "No NetCDF trajectory reporter available: neither "
                    "openmm.app.NetCDFReporter nor mdtraj is importable. The MD "
                    "trajectory cannot be written for MMPBSA.py. Install mdtraj: "
                    "conda install -n <env> -c conda-forge mdtraj")
        else:
            errors.append(
                "OpenMM is not importable in the probe env -- MD cannot run. "
                f"Detail: {payload.get('error','unknown')}")
            if cp.stderr.strip():
                info["probe_stderr"] = cp.stderr.strip()[-800:]
    except Exception as e:  # noqa: BLE001
        warnings.append(f"Could not run the OpenMM probe ({type(e).__name__}: "
                        f"{e}); proceeding without it.")

    return {"ok": not errors, "errors": errors,
            "warnings": warnings, "info": info}


def _summarize_unit_failures(data) -> dict:
    """Mine an OpenFE result's unit_results for ProtocolUnitResult / *Failure.

    Returns {"n_units","n_ok","n_failed","failures":[...], "headline": str}.
    Each failure entry: {"unit","exc_type","exc_msg","origin","traceback"}.
    `origin` is the deepest stack frame that lives in site-packages of the
    failing library (so we point at the real culprit, not the generic wrapper).
    """
    ur = data.get("unit_results", {})
    items = list(ur.values()) if isinstance(ur, dict) else list(ur or [])
    n_units = len(items)
    failures, n_ok = [], 0
    for u in items:
        if not isinstance(u, dict):
            continue
        qn = u.get("__qualname__", "")
        exc = u.get("exception")
        tb = u.get("traceback") or ""
        if not exc and "Failure" not in qn:
            n_ok += 1
            continue
        # exception is typically ["ExcType", [args...]] or a string
        exc_type, exc_msg = "", ""
        if isinstance(exc, (list, tuple)) and exc:
            exc_type = str(exc[0])
            rest = exc[1] if len(exc) > 1 else ""
            if isinstance(rest, (list, tuple)):
                exc_msg = ", ".join(str(x) for x in rest)
            else:
                exc_msg = str(rest)
        elif isinstance(exc, str):
            exc_type = exc
        # origin: deepest 'File ".../site-packages/<pkg>/..." line N, in func'
        origin = ""
        if tb:
            import re as _re
            frames = _re.findall(
                r'File "([^"]+)", line (\d+), in (\S+)', tb)
            # prefer frames inside site-packages but not gufe's generic wrapper
            picks = [f for f in frames if "site-packages" in f[0]
                     and "protocolunit.py" not in f[0]]
            chosen = (picks[-1] if picks else (frames[-1] if frames else None))
            if chosen:
                fn = chosen[0].split("site-packages/")[-1]
                origin = f"{fn}:{chosen[1]} in {chosen[2]}()"
        failures.append({
            "unit": u.get("name", qn or "unit"),
            "exc_type": exc_type, "exc_msg": exc_msg,
            "origin": origin, "traceback": tb})
    n_failed = len(failures)
    headline = ""
    if n_failed:
        f0 = failures[0]
        loc = f" in {f0['origin']}" if f0["origin"] else ""
        # de-duplicate: how many failures share this exc_type
        same = sum(1 for f in failures if f["exc_type"] == f0["exc_type"])
        mult = f" x{same}" if same > 1 else ""
        msg = (f": {f0['exc_msg']}" if f0["exc_msg"] else "")
        headline = f"{f0['exc_type'] or 'failure'}{mult}{msg}{loc}"
    return {"n_units": n_units, "n_ok": n_ok, "n_failed": n_failed,
            "failures": failures, "headline": headline}


def validate_openfe_result(result_json) -> dict:
    """Inspect an OpenFE *_result.json and classify it honestly.

    Returns {"status": "completed"|"failed"|"unreadable",
             "estimate": float|None, "reason": str, "size_bytes": int,
             "failure": {...}|None}.

    Key fixes:
      * #2 -- a JSON with estimate=null is NOT "completed".
      * failure mining -- if the run failed, dig the REAL exception/traceback
        out of unit_results (ProtocolUnitFailure) instead of just saying
        "no estimate". `failure` carries the per-unit details + a one-line
        `reason` headline (e.g. the Numba JIT FileNotFoundError).
    """
    import json
    from pathlib import Path
    p = Path(result_json)
    info = {"status": "failed", "estimate": None, "reason": "",
            "size_bytes": p.stat().st_size if p.exists() else 0,
            "failure": None}
    if not p.exists():
        info["status"] = "unreadable"; info["reason"] = "missing"; return info
    try:
        d = json.loads(p.read_text())
    except Exception as e:  # noqa: BLE001
        info["status"] = "unreadable"; info["reason"] = f"unparseable JSON: {e}"
        return info

    est = d.get("estimate")
    if isinstance(est, dict):
        est = est.get("magnitude")

    fsum = _summarize_unit_failures(d)

    # Top-level explicit failure markers (some protocols write these).
    top_exc = d.get("exception") or d.get("traceback") or d.get("error")

    if est is not None:
        info["status"] = "completed"; info["estimate"] = est
        info["reason"] = "estimate present"
        if fsum["n_failed"]:
            # estimate present but some replicates/units failed -> warn
            info["reason"] += (f" (warning: {fsum['n_failed']} unit(s) failed: "
                               f"{fsum['headline']})")
            info["failure"] = fsum
        return info

    # No estimate -> failed. Surface the real reason.
    info["status"] = "failed"
    if fsum["n_failed"]:
        info["failure"] = fsum
        info["reason"] = (f"{fsum['n_failed']}/{fsum['n_units']} unit(s) failed "
                          f"-- {fsum['headline']}")
    elif top_exc:
        info["reason"] = "result carries a top-level exception/traceback"
        info["failure"] = {"n_units": fsum["n_units"], "n_ok": fsum["n_ok"],
                           "n_failed": 0, "failures": [],
                           "headline": "top-level exception",
                           "top_traceback": (top_exc if isinstance(top_exc, str)
                                             else str(top_exc))}
    else:
        info["reason"] = "no estimate and no protocol payload"
    return info


def scan_results_dir(root) -> dict:
    """Disk-scan a work_dir for finished jobs, independent of session_state (#1).

    Detects BOTH workflows:
      * OpenFE FEP edges : <root>/**/*_result.json  (validated)
      * MM-GBSA screens  : <root>/**/lig_*/mmgbsa/FINAL_RESULTS_MMPBSA.dat

    Returns {"rows": [ {kind,name,status,detail,path,estimate} ],
             "counts": {status: n}}.
    """
    import re
    from pathlib import Path
    root = Path(root).expanduser()
    rows = []
    if not root.exists():
        return {"rows": rows, "counts": {}, "error": f"not found: {root}"}

    # --- OpenFE result JSONs ---
    for rj in sorted(root.rglob("*_result.json")):
        v = validate_openfe_result(rj)
        est = v["estimate"]
        if isinstance(est, (int, float)):
            detail = f"DG/DDG={est:.2f}  ({v['reason']}; {v['size_bytes']:,} B)"
        else:
            detail = f"{v['reason']}; {v['size_bytes']:,} B"
        rows.append({"kind": "OpenFE-result", "name": rj.stem,
                     "status": v["status"], "detail": detail,
                     "path": str(rj), "estimate": est,
                     "failure": v.get("failure")})

    # --- OpenFE edges that have a _work dir but NO result yet (running/pending) ---
    seen_results = {Path(r["path"]).name.replace("_result.json", "")
                    for r in rows}
    for wd in sorted(root.rglob("*_work")):
        stem = wd.name[:-5]  # strip '_work'
        if stem in seen_results:
            continue
        rows.append({"kind": "OpenFE-edge", "name": stem,
                     "status": "running",
                     "detail": "work dir present, no result.json yet",
                     "path": str(wd), "estimate": None})

    # --- MM-GBSA per-ligand results ---
    for dat in sorted(root.rglob("lig_*/mmgbsa/FINAL_RESULTS_MMPBSA.dat")):
        lig = dat.parent.parent.name  # lig_<name>
        dg = None
        try:
            m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)", dat.read_text())
            if m:
                dg = float(m.group(1))
        except Exception:  # noqa: BLE001
            pass
        detail = (f"DG_bind={dg:.2f} kcal/mol" if dg is not None
                  else "result file present (unparsed)")
        rows.append({"kind": "MM-GBSA", "name": lig, "status": "completed",
                     "detail": detail, "path": str(dat), "estimate": dg})

    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    return {"rows": rows, "counts": counts, "error": None}
