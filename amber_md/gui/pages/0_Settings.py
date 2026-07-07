"""
0_Settings.py  --  Site / Environment configuration.

Edit everything that changes when this workflow moves to a new machine or
cluster in ONE place: scheduler (LSF now; SLURM planned), GPU/CPU resources,
Amber and OpenFE locations, and default paths. Values are persisted to a site
config file (YAML/JSON) and read by the other pages as their field defaults.

Resolution order: $BFEP_SITE_CONFIG -> ~/.bfep/site_config.* -> <repo>/site_config.* -> built-in defaults.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
# -*- coding: utf-8 -*-
from __future__ import annotations

import shutil
from dataclasses import asdict

import streamlit as st

st.set_page_config(page_title="Settings", layout="wide", page_icon="gear")

from amber_md import site_config as sc

st.title("Settings - Site & Environment")
st.caption("Adapt this workflow to your environment here. Other pages read these "
           "values as their defaults. Nothing here changes the science - only "
           "where jobs run and how tools are located.")

cfg = sc.load()
_path = sc.resolve_config_path()
if _path:
    st.info(f"Loaded site config from: `{_path}`")
else:
    st.warning("No site config file found yet - showing built-in defaults. "
               "Save below to create one.")

# ---- Scheduler ------------------------------------------------------------
st.header("1 - Scheduler")
c1, c2, c3 = st.columns(3)
sched_type = c1.selectbox(
    "Scheduler type", ["lsf", "slurm", "local"],
    index=["lsf", "slurm", "local"].index(cfg.scheduler.type)
    if cfg.scheduler.type in ("lsf", "slurm", "local") else 0,
    help="LSF is fully supported. SLURM is a stub (not yet validated). "
         "local runs without a batch scheduler.")
if sched_type == "slurm":
    st.warning("SLURM is not yet validated end-to-end. Submissions will refuse "
               "until amber_md/schedulers/slurm.py is implemented + tested.")
submit_cmd = c2.text_input("Submit command", cfg.scheduler.submit_cmd)
query_cmd = c3.text_input("Query command", cfg.scheduler.query_cmd)
c1, c2, c3 = st.columns(3)
cancel_cmd = c1.text_input("Cancel command", cfg.scheduler.cancel_cmd)
gpu_queue = c2.text_input("GPU queue / partition", cfg.scheduler.gpu_queue)
cpu_queue = c3.text_input("CPU queue / partition", cfg.scheduler.cpu_queue)
c1, c2, c3 = st.columns(3)
project = c1.text_input("Project / account", cfg.scheduler.project)
walltime = c2.text_input("Default walltime", cfg.scheduler.walltime)
max_conc = c3.number_input("Max concurrent jobs", 1, 512,
                           int(cfg.scheduler.max_concurrent))

# ---- Resources ------------------------------------------------------------
st.header("2 - GPU / CPU resources")
c1, c2, c3 = st.columns(3)
n_gpu = c1.number_input("GPUs per job", 1, 16, int(cfg.resources.n_gpu))
gpu_mode = c2.selectbox("GPU request mode", ["rusage", "slots"],
                        index=0 if cfg.resources.gpu_request_mode == "rusage" else 1)
gpu_cpu_cores = c3.number_input("CPU cores alongside GPU", 1, 128,
                                int(cfg.resources.gpu_cpu_cores))
c1, c2, c3 = st.columns(3)
fep_gpu_cores = c1.number_input("FEP window cores", 1, 128,
                                int(cfg.resources.fep_gpu_cores))
fep_mem_mb = c2.number_input("FEP mem per window (MB)", 512, 262144,
                             int(cfg.resources.fep_mem_mb), 512)
n_cpu = c3.number_input("MM-GBSA MPI ranks", 1, 256, int(cfg.resources.n_cpu))
c1, c2 = st.columns(2)
cap_gpu = c1.number_input("CPU-settle cap (GPU nodes)", 1, 256,
                          int(cfg.resources.cpu_settle_cap_gpu))
cap_cpu = c2.number_input("CPU-settle cap (CPU queue)", 1, 256,
                          int(cfg.resources.cpu_settle_cap_cpu))
avoid = st.text_input("Avoid hosts (comma-separated)",
                      ", ".join(cfg.resources.avoid_hosts),
                      help="Nodes with known driver/ECC issues, excluded from "
                           "every GPU job.")

# ---- Amber ----------------------------------------------------------------
st.header("3 - Amber")
c1, c2 = st.columns(2)
amberhome = c1.text_input("AMBERHOME (optional)", cfg.amber.amberhome,
                          help="Leave blank to rely on the environment / modules.")
pmemd = c2.text_input("pmemd.cuda (binary or path)", cfg.amber.pmemd_cuda)
modules = st.text_input("module load (space-separated)",
                        " ".join(cfg.amber.module_load))
charge_method = st.selectbox("Default charge method", ["bcc", "gas", "resp"],
                             index=["bcc", "gas", "resp"].index(cfg.amber.charge_method)
                             if cfg.amber.charge_method in ("bcc", "gas", "resp") else 0)

# ---- OpenFE ---------------------------------------------------------------
st.header("4 - OpenFE / OpenMM")
c1, c2, c3 = st.columns(3)
ofe_env = c1.text_input("conda env", cfg.openfe.conda_env)
ofe_py = c2.text_input("python bin (optional)", cfg.openfe.python_bin,
                       help="Explicit interpreter; overrides conda env if set.")
smff = c3.text_input("Small-molecule FF", cfg.openfe.small_molecule_ff)

# ---- Paths ----------------------------------------------------------------
st.header("5 - Paths")
c1, c2, c3 = st.columns(3)
run_dir = c1.text_input("Default run dir", cfg.paths.default_run_dir)
scratch = c2.text_input("Scratch dir (optional)", cfg.paths.scratch_dir)
venv = c3.text_input("venv activate script", cfg.paths.venv_activate)

# ---- Save -----------------------------------------------------------------
st.divider()
col_a, col_b = st.columns([1, 2])
save_to = col_b.text_input(
    "Save to", str(sc.default_save_path()),
    help="Per-user by default (~/.bfep/). Point at the repo's site_config.yaml "
         "to share a site-wide default.")
if col_a.button("Save settings", type="primary", width="stretch"):
    new = sc.SiteConfig.from_dict({
        "scheduler": {"type": sched_type, "submit_cmd": submit_cmd,
                      "query_cmd": query_cmd, "cancel_cmd": cancel_cmd,
                      "gpu_queue": gpu_queue, "cpu_queue": cpu_queue,
                      "project": project, "walltime": walltime,
                      "max_concurrent": int(max_conc)},
        "resources": {"n_gpu": int(n_gpu), "gpu_request_mode": gpu_mode,
                      "gpu_cpu_cores": int(gpu_cpu_cores),
                      "fep_gpu_cores": int(fep_gpu_cores),
                      "fep_mem_mb": int(fep_mem_mb), "n_cpu": int(n_cpu),
                      "cpu_settle_cap_gpu": int(cap_gpu),
                      "cpu_settle_cap_cpu": int(cap_cpu),
                      "avoid_hosts": [h.strip() for h in avoid.split(",") if h.strip()]},
        "amber": {"amberhome": amberhome, "pmemd_cuda": pmemd,
                  "module_load": modules.split(), "charge_method": charge_method},
        "openfe": {"conda_env": ofe_env, "python_bin": ofe_py,
                   "small_molecule_ff": smff},
        "paths": {"default_run_dir": run_dir, "scratch_dir": scratch,
                  "venv_activate": venv},
    })
    try:
        written = sc.save(new, save_to)
        sc.refresh()
        st.success(f"Saved site config to `{written}`. Reload the other pages "
                   "to pick up the new defaults.")
    except Exception as e:  # noqa: BLE001
        st.error(f"Save failed: {e}")

# ---- Validate -------------------------------------------------------------
st.divider()
st.header("6 - Test environment")
if st.button("Run checks"):
    rows = []
    def _check(label, ok, detail=""):
        rows.append({"check": label, "status": "OK" if ok else "MISSING",
                     "detail": detail})
    _check(f"scheduler submit ({submit_cmd})", shutil.which(submit_cmd) is not None,
           shutil.which(submit_cmd) or "not on PATH")
    _check(f"pmemd ({pmemd})", shutil.which(pmemd) is not None or "/" in pmemd,
           shutil.which(pmemd) or pmemd)
    import os as _os
    _check("AMBERHOME", bool(amberhome) or bool(_os.environ.get("AMBERHOME")),
           amberhome or _os.environ.get("AMBERHOME", ""))
    _check(f"conda ({ofe_env})", shutil.which("conda") is not None,
           shutil.which("conda") or "conda not on PATH")
    import pandas as pd
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    st.caption("This is a lightweight PATH check. For a full diagnostic run "
               "`python check_env.py` in the target environment.")
