"""
3_MMGBSA_Screen.py  --  MM-GBSA multi-ligand screen: submit / resume / aggregate.

A METHOD-named screen page (not engine-named). It manages an MM-GBSA batch run
for BOTH engines:

  * Amber  : python -m amber_md.batch  (pmemd MD -> MMPBSA.py), with resume.
  * OpenMM : per-ligand fan-out is launched from Setup & Launch; this page still
             handles its AGGREGATION (both engines write
             lig_*/mmgbsa/FINAL_RESULTS_MMPBSA.dat, so one aggregator serves both).

Scope is deliberately limited to the verbs that **Results - Compare** does NOT
provide:
  * Submit  (Amber batch screen; OpenMM screens are submitted from page 0)
  * Resume / status      (python -m amber_md.batch_resume)
  * Aggregate            (python -m amber_md.batch_aggregate -> INDEX.html)

RANKING / CHARTING / EXPERIMENTAL RMSE is NOT duplicated here: the
**Results - Compare** page is the universal ranker -- it already auto-detects and
ranks MM-GBSA, ABFE, and RBFE results across a campaign. Use this page to run a
screen; use Results - Compare to rank it.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="MM-GBSA Screen", layout="wide", page_icon="rocket")

from amber_md.gui.common import dir_picker, file_picker, spawn_detached
from amber_md import site_config as _sc
_SITE = _sc.get()

st.title("MM-GBSA Screen - Submit, Resume, Aggregate")
st.caption("Manage a multi-ligand MM-GBSA screen (Amber or OpenMM). "
           "For ranking, charts, and experimental RMSE, use **Results - Compare** "
           "(it ranks MM-GBSA, ABFE, and RBFE results across the campaign).")

exp = st.session_state.get("experiment", {})
default_wd = exp.get("work_dir", "~/Run_dir/screen")

with st.sidebar:
    st.header("Batch directory")
    bdir = dir_picker("Batch dir (--batch-dir)", "mms_bdir",
                      default_dir=default_wd)
    st.divider()
    st.header("Engine")
    engine = st.radio(
        "MD engine for SUBMIT", ["Amber", "OpenMM"], 0, key="mms_engine",
        help="Amber screens submit here. OpenMM screens are submitted from "
             "Setup & Launch (per-ligand fan-out); use this page only to "
             "Aggregate them.")
    st.divider()
    st.header("Submit settings (Amber)")
    queue = st.text_input("GPU queue", _SITE.scheduler.gpu_queue, key="mms_queue")
    walltime = st.text_input("Walltime", _SITE.scheduler.walltime, key="mms_wall")
    project = st.text_input("Project / account (optional)", _SITE.scheduler.project, key="mms_proj")
    maxc = st.number_input("Max concurrent jobs", 1, 512, int(_SITE.scheduler.max_concurrent), key="mms_maxc")
    prod_ns = st.number_input("Production (ns)", 0.1, 200.0, 50.0, 0.5,
                              key="mms_prod")
    equil_ns = st.number_input("Equilibration (ns)", 0.1, 20.0, 1.0, 0.1,
                               key="mms_equil")
    charge_method = st.selectbox("Charge method", ["bcc", "gas", "resp"], ["bcc","gas","resp"].index(_SITE.amber.charge_method) if _SITE.amber.charge_method in ("bcc","gas","resp") else 0,
                                 key="mms_chg")
    salt = st.number_input("Salt conc (M)", 0.0, 1.0, 0.15, 0.01, key="mms_salt")
    no_gbsa = st.checkbox("Skip MM-GBSA (MD only)", value=False, key="mms_nogbsa")
    decomp = st.checkbox("Per-residue decomposition", value=False, key="mms_dec")
    decomp_mask = st.text_input("Decomp residue mask (e.g. ':300-450')", "",
                                key="mms_decmask", disabled=not decomp)
    dry = st.checkbox("Dry run (scaffold only, no bsub)", value=False,
                      key="mms_dry")

bdir_p = Path(bdir).expanduser() if bdir else None
PY = sys.executable

# ---- Inputs (only needed for Amber submit) ----
st.subheader("1 - Inputs (for Amber submit)")
ci1, ci2 = st.columns(2)
with ci1:
    protein = file_picker("Protein (--protein-file)", "mms_protein",
                          {".pdb", ".mol2", ".cif"}, default_dir="~/Run_dir")
with ci2:
    ligands = file_picker("Ligands (multi-record SDF/MOL2 or dir)",
                          "mms_ligands", {".sdf", ".mol2", ".mol"},
                          default_dir="~/Run_dir")

# ---- Actions ----
st.subheader("2 - Actions")
a1, a2, a3 = st.columns(3)
amber_submit_ok = (engine == "Amber" and protein and ligands and bdir_p)
submit = a1.button("Submit Amber screen", type="primary",
                   disabled=not amber_submit_ok, key="mms_submit")
resume = a2.button("Resume / status (Amber)", disabled=not bdir_p,
                   key="mms_resume")
aggregate = a3.button("Aggregate + rank (both engines)", disabled=not bdir_p,
                      key="mms_agg")

if engine == "OpenMM":
    st.info("OpenMM MM-GBSA screens are submitted from **Setup & Launch** "
            "(choose MM-GBSA + OpenMM/OpenFE; it fans out one job per ligand). "
            "Come back here to **Aggregate + rank** once they finish -- the "
            "aggregator reads `lig_*/mmgbsa/FINAL_RESULTS_MMPBSA.dat`, which "
            "the OpenMM fan-out writes just like the Amber path.")

if submit:
    bdir_p.mkdir(parents=True, exist_ok=True)
    cmd = [
        PY, "-m", "amber_md.batch",
        "--protein-file", str(Path(protein)),
        "--ligands", str(Path(ligands)),
        "--batch-dir", str(bdir_p),
        "--queue", str(queue),
        "--walltime", str(walltime),
        "--prod-ns", str(prod_ns),
        "--equil-ns", str(equil_ns),
        "--charge-method", str(charge_method),
        "--salt", str(salt),
        "--max-concurrent", str(int(maxc)),
    ]
    if project.strip():
        cmd += ["--project", project.strip()]
    if no_gbsa:
        cmd.append("--no-gbsa")
    if decomp:
        cmd.append("--decomp")
        if decomp_mask.strip():
            cmd += ["--decomp-residues", decomp_mask.strip()]
    if dry:
        cmd.append("--dry-run")
    log = bdir_p / "mmgbsa_screen_submit.log"
    st.code(" ".join(cmd), language="bash")
    try:
        pid = spawn_detached(cmd, log, cwd=str(bdir_p))
        st.success(f"Amber screen launched (pid {pid}). Log: {log}")
        st.session_state["experiment"] = {**exp, "work_dir": str(bdir_p),
                                          "method": "MM-GBSA", "engine": "Amber"}
    except Exception as e:  # noqa: BLE001
        st.error(f"Launch failed: {e}")

if resume:
    cmd = [PY, "-m", "amber_md.batch_resume", str(bdir_p)]
    log = bdir_p / "mmgbsa_screen_resume.log"
    st.code(" ".join(cmd), language="bash")
    try:
        pid = spawn_detached(cmd, log, cwd=str(bdir_p))
        st.success(f"Resume/status launched (pid {pid}). Log: {log}")
    except Exception as e:  # noqa: BLE001
        st.error(f"Launch failed: {e}")

if aggregate:
    cmd = [PY, "-m", "amber_md.batch_aggregate", str(bdir_p)]
    log = bdir_p / "mmgbsa_screen_aggregate.log"
    st.code(" ".join(cmd), language="bash")
    try:
        pid = spawn_detached(cmd, log, cwd=str(bdir_p))
        st.success(f"Aggregate launched (pid {pid}). Log: {log}. "
                   "Open INDEX.html when it finishes -- or use "
                   "**Results - Compare** for an interactive ranked view.")
    except Exception as e:  # noqa: BLE001
        st.error(f"Launch failed: {e}")
