"""
1_Setup_and_Launch.py  --  Unified Setup & Launch wizard (v2.5.0, Phase 1).

A single page that replaces the method/engine-specific setup pages with three
selectors (Method x Engine x Scope) plus shared inputs. It auto-detects ligand
count, enforces a Method x Engine compatibility matrix, lets the user review a
method-specific parameter panel, and LAUNCHES by calling the existing, unchanged
engine functions:

  * RBFE / OpenFE  -> openfe_common.plan_rbfe_cmd  + spawn_detached
  * ABFE / OpenFE  -> python -m amber_md.abfe_openfe_plan (via OpenFE env python)
  * MM-GBSA / Amber -> amber_md.run_amber (CLI) via spawn_detached

NOTHING in the science/engine layer changes. This page only orchestrates.

The selected experiment is stored in st.session_state["experiment"] so the
Job Monitor and Results pages (Phases 2-3) can read protein/ligand/engine/method
without re-prompting.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from datetime import datetime as _dt

import shlex
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="Setup & Launch", layout="wide", page_icon="rocket")

from amber_md.gui.common import (
    file_picker, spawn_detached, sanitized_openfe_env, run_amber_py, get_lsf_jobs,
)
from amber_md.gui import amber_config
from amber_md.gui.fep_common import load_sdf_records
from amber_md.gui.openfe_common import (
    OpenFESettings, write_network_yaml, plan_rbfe_cmd,
    preflight_openfe_charges, preflight_openmm_mmgbsa,
)
from amber_md import site_config as _sc
_SITE = _sc.get()

# ----------------------------------------------------------------------------
# Static config: the Method x Engine compatibility matrix.
# value = (enabled, note). Disabled combos are shown greyed with the note.
# ----------------------------------------------------------------------------
METHODS = ["MM-GBSA", "ABFE", "RBFE"]
ENGINES = ["Amber", "OpenMM / OpenFE"]

# NOTE: Amber ABFE and Amber RBFE were removed from the supported GUI surface
# (they were never validated end-to-end). The underlying engine code
# (amber_md.fep_driver, rbfe_map, abfe_* modules) is retained in the package for
# programmatic / CLI use, but these two combinations are intentionally disabled
# here so the GUI exposes only supported, functional workflows:
#     * MM-GBSA / Amber            (native AmberTools pipeline)
#     * MM-GBSA / OpenMM / OpenFE  (OpenMM MD + AmberTools scoring)
#     * ABFE    / OpenMM / OpenFE  (OpenFE AbsoluteBindingProtocol)
#     * RBFE    / OpenMM / OpenFE  (OpenFE relative binding network)
COMPAT = {
    ("MM-GBSA", "Amber"):           (True,  "Native AmberTools MMPBSA.py pipeline."),
    ("MM-GBSA", "OpenMM / OpenFE"): (True,  "OpenMM MD + AmberTools MMPBSA.py (experimental)."),
    ("ABFE",    "Amber"):           (False, "Amber ABFE is not exposed in the GUI. "
                                            "Use ABFE / OpenMM / OpenFE."),
    ("ABFE",    "OpenMM / OpenFE"): (True,  "OpenFE AbsoluteBindingProtocol."),
    ("RBFE",    "Amber"):           (False, "Amber RBFE is not exposed in the GUI. "
                                            "Use RBFE / OpenMM / OpenFE."),
    ("RBFE",    "OpenMM / OpenFE"): (True,  "OpenFE relative binding network."),
}

SCOPES = ["Auto-detect", "Single ligand", "Multi-ligand"]


@dataclass
class Experiment:
    method: str
    engine: str
    scope: str
    protein: str | None
    ligands: str | None
    work_dir: str
    n_ligands: int
    effective_scope: str           # resolved single/multi after auto-detect
    params: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def _segmented(label: str, options: list[str], key: str, default: str | None = None):
    """Segmented control with graceful fallback to radio on older Streamlit."""
    idx = options.index(default) if default in options else 0
    if hasattr(st, "segmented_control"):
        val = st.segmented_control(label, options, default=options[idx], key=key)
        return val or options[idx]
    return st.radio(label, options, index=idx, key=key, horizontal=True)


def _resolve_scope(scope: str, n: int) -> str:
    if scope == "Single ligand":
        return "single"
    if scope == "Multi-ligand":
        return "multi"
    # auto
    return "multi" if n >= 2 else "single"


def _compat(method: str, engine: str) -> tuple[bool, str]:
    return COMPAT.get((method, engine), (False, "Unsupported combination."))


# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------
st.title("Setup & Launch")
st.caption("Pick a method, an engine, and your inputs. The page adapts and "
           "launches the right workflow -- one place for everything.")

with st.expander("Need help choosing? (method / engine / scope)", expanded=False):
    st.markdown("""
**Method** - what kind of calculation:
- **MM-GBSA** - fast, approximate ranking. Best for triaging many ligands.
- **ABFE** - rigorous *absolute* binding free energy per ligand. Slower, accurate.
- **RBFE** - rigorous *relative* free energy across a series of similar ligands
  (needs >= 2 ligands). Slower, accurate.

**Engine** - the software that runs it:
- **Amber** - native AmberTools MM-GBSA.
- **OpenMM / OpenFE** - OpenMM MD; the only engine for ABFE / RBFE here.

**Scope** - how many ligands:
- **Auto-detect** (recommended) - reads the count from your ligand file.
- **Single / Multi-ligand** - force one or a batch.

*A green banner below confirms a supported combination; red means pick a
different engine (it suggests one).*
""")

# ============================ 1. Selectors ============================
st.subheader("1 - What do you want to run?")
c1, c2, c3 = st.columns(3)
with c1:
    method = _segmented("Method", METHODS, "wiz_method", "ABFE")
with c2:
    engine = _segmented("Engine", ENGINES, "wiz_engine", "OpenMM / OpenFE")
with c3:
    scope_choice = _segmented("Scope", SCOPES, "wiz_scope", "Auto-detect")

ok, note = _compat(method, engine)
if ok:
    st.success(f"**{method} / {engine}** - {note}")
else:
    st.error(f"**{method} / {engine}** is not available. {note}")
    # Suggest the nearest valid engine for this method
    alts = [e for e in ENGINES if _compat(method, e)[0]]
    if alts:
        st.info(f"Tip: for **{method}**, use **{alts[0]}**.")

# ============================ 2. Inputs ============================
st.subheader("2 - Inputs")
# final42 FIX: `params` is written by the ligand-residue-name widget below
# (well before the "3 - Parameters" section), but it was only initialized in
# section 3 -> NameError: name 'params' is not defined on first render.
# Initialize it here so early writers (ligand_resname) are safe.
params: dict = {}
ci1, ci2 = st.columns(2)
with ci1:
    protein = file_picker("Protein (PDB / MOL2 / CIF)", "wiz_protein",
                          {".pdb", ".mol2", ".cif", ".pdbx"},
                          default_dir="~/Run_dir")
with ci2:
    ligands = file_picker("Ligands (SDF / MOL2 -- 1 or many records)",
                          "wiz_ligands", {".sdf", ".mol", ".mol2"},
                          default_dir="~/Run_dir")

    # Ligand residue name: auto-detected from a .mol2 (e.g. 'UNK'); editable.
    # Empty => the engine auto-detects (mol2) or falls back to 'LIG' (sdf).
    # Drives ABFE alchemical masks (timask/scmask/crgmask); shown for RBFE /
    # MM-GBSA as the residue tleap/antechamber will assign.
    _auto_rn, _rn_src = ("", "")
    if ligands:
        try:
            from amber_md.utils import detect_ligand_resname as _detect_rn
            _auto_rn, _rn_src = _detect_rn(ligands)
        except Exception:
            _auto_rn, _rn_src = ("", "")
    _rn_help = ("Ligand residue name. Leave blank to auto-detect from a .mol2 "
                "(SDF -> assigned by tleap as 'LIG'). For ABFE this sets the "
                "alchemical masks; for RBFE/MM-GBSA it is informational.")
    params["ligand_resname"] = st.text_input(
        "Ligand residue name (blank = auto)", value="", key="wiz_resname",
        help=_rn_help)
    if ligands and _rn_src == "mol2":
        _eff = params["ligand_resname"].strip() or _auto_rn
        st.caption(f"Detected residue name in mol2: **{_auto_rn}**  -> "
                   f"effective: **{_eff}**"
                   + ("" if params['ligand_resname'].strip()
                      else "  (auto-detected; type above to override)"))
    elif ligands:
        st.caption("No residue name in this file (SDF) -- tleap will assign "
                   "'LIG' unless you set one above.")

work_dir = st.text_input("Work directory", "~/Run_dir/run_v250",
                         key="wiz_wd")

# ---- Execution target (applies to all methods) ----
_et1, _et2, _et3 = st.columns([1.2, 1, 1])
exec_target = _et1.selectbox(
    "Execution target",
    ["GPU queue (LSF)", "Local host"],
    index=0, key="wiz_exec_target",
    help="GPU queue submits to the cluster GPU queue via the configured "
         "scheduler (Settings page) -- the default. Local host runs on THIS "
         "machine's GPU. Applies to MM-GBSA, RBFE and ABFE.")
_is_local = (exec_target == "Local host")
gpu_queue = _et2.text_input("GPU queue", _SITE.scheduler.gpu_queue,
                            key="wiz_gpu_queue", disabled=_is_local)
gpu_walltime = _et3.text_input("Walltime (HH:MM)", _SITE.scheduler.walltime,
                               key="wiz_gpu_walltime", disabled=_is_local)
# ---- Optional submission throttle (login-node, survives page close) ----
_tt1, _tt2 = st.columns([1.2, 1])
throttle_on = _tt1.checkbox(
    "Throttle submissions", value=False, key="wiz_throttle_on",
    disabled=_is_local,
    help="OFF (default): all jobs are submitted at once. ON: a small detached "
         "helper on the login node keeps no more than N of THIS batch's jobs in "
         "the queue (PEND+RUN) at a time, submitting the next as slots free. It "
         "runs under nohup, so it survives closing this page. Useful for large "
         "RBFE/ABFE campaigns so you don't flood the shared queue.")
throttle_n = _tt2.number_input(
    "Max jobs in queue", min_value=1, max_value=200, value=8, step=1,
    key="wiz_throttle_n", disabled=(_is_local or not throttle_on))

# ---- auto-detect ligand count ----
n_lig = 0
lig_records = []
if ligands:
    lig_records = load_sdf_records(Path(ligands))
    n_lig = len(lig_records)

eff_scope = _resolve_scope(scope_choice, n_lig)

# input status line
ic1, ic2, ic3 = st.columns(3)
ic1.metric("Protein", "OK" if protein else "—",
           Path(protein).name if protein else "not set")
ic2.metric("Ligands detected", n_lig if ligands else "—",
           Path(ligands).name if ligands else "not set")
ic3.metric("Resolved scope", eff_scope.upper() if ligands else "—",
           f"from '{scope_choice}'")

if ligands and n_lig == 0:
    st.warning("Could not parse any ligand records. Check the file / RDKit.")
if ligands and lig_records:
    with st.expander(f"Ligand records ({n_lig})", expanded=False):
        st.dataframe(
            [{"#": r.index, "name": r.name, "heavy_atoms": r.n_heavy,
              "smiles": (r.smiles or "")[:60]} for r in lig_records],
            hide_index=True, width="stretch")

# scope/method sanity hints
if method == "RBFE" and eff_scope == "single":
    st.warning("RBFE needs >=2 ligands to form a perturbation network. "
               "Provide a multi-record SDF or switch method.")
if method in ("MM-GBSA",) and eff_scope == "multi":
    st.info(f"Multi-ligand MM-GBSA -> batch ranking of {n_lig} ligands.")
if method == "RBFE" and eff_scope == "multi":
    st.info(f"Multi-ligand RBFE -> minimal spanning network over {n_lig} ligands.")

# ============================ 3. Parameters ============================
st.subheader("3 - Parameters")
# final42 FIX: do NOT re-bind params here (it was initialized in section 2 and
# may already hold ligand_resname); just ensure it exists.
try:
    params  # noqa: F821  (defined in section 2)
except NameError:
    params: dict = {}

if method == "MM-GBSA":
    with st.expander("MM-GBSA settings", expanded=True):
        p1, p2, p3 = st.columns(3)
        params["igb"] = p1.selectbox("GB model (igb)", [8, 5, 2, 1], 0,
                                     key="wiz_igb")
        params["salt_conc"] = p2.number_input("Salt conc (M)", 0.0, 1.0, 0.15,
                                              0.01, key="wiz_salt")
        params["stride"] = p3.number_input("Frame stride", 1, 100, 1,
                                           key="wiz_stride")
        params["decomposition"] = st.checkbox(
            "Per-residue decomposition (FINAL_DECOMP_MMPBSA.dat)",
            value=False, key="wiz_decomp")
        if params["decomposition"]:
            params["decomp_residues"] = st.text_input(
                "Decomp residue mask (Amber, e.g. ':1-50' or blank=all)",
                "", key="wiz_decompres")
        # final43: expose MM-GBSA MD simulation time for BOTH engines (was
        # hard-coded and unreachable from the GUI). final51: unified to
        # 10 ns for BOTH engines so Amber and OpenMM scores are comparable.
        _is_amber_mm = (engine == "Amber")
        _def_prod = 10.0  # final51: unified 10 ns for BOTH engines
        st.markdown("**MD simulation time** (applies to the selected engine)")
        t1, t2 = st.columns(2)
        params["prod_ns"] = t1.number_input(
            "Production MD (ns)", 0.1, 200.0, _def_prod, 0.1,
            key="wiz_mmgbsa_prodns",
            help="MD production length sampled for MM-GBSA frames. Default "
                 "10 ns for both Amber and OpenMM engines. Longer = better "
                 "averaging, more GPU time.")
        params["equil_ns"] = t2.number_input(
            "Equilibration MD (ns)", 0.0, 20.0, 1.0, 0.1, key="wiz_mmgbsa_eqns",
            help="Equilibration before production (default 1 ns); not sampled "
                 "for MM-GBSA.")
        st.markdown("**OpenMM engine options** (used only when Engine = OpenMM / OpenFE)")
        o1, o2, o3 = st.columns(3)
        params["platform"] = o1.selectbox(
            "OpenMM platform", ["CUDA", "OpenCL", "CPU"], 0, key="wiz_omm_plat")
        params["buffer_A"] = o2.number_input(
            "Solvent box buffer (A)", 6.0, 16.0, 12.0, 0.5, key="wiz_omm_buf")  # v2.5.71: default 10->12 (match CLI/config)
        params["mm_charge_method"] = o3.selectbox(
            "Ligand charge method", ["bcc", "gas"], 0, key="wiz_omm_chg")
        params["ligand_charge"] = st.number_input(
            "Ligand net charge", -5, 5, 0, 1, key="wiz_omm_lchg")
        fr1, fr2, fr3 = st.columns(3)
        params["mmgbsa_start_frame"] = fr1.number_input(
            "Start frame", 1, 100000, 1, key="wiz_startf")
        params["mmgbsa_end_frame"] = fr2.number_input(
            "End frame (0 = last)", 0, 100000, 0, key="wiz_endf")
        params["mmgbsa_stride"] = fr3.number_input(
            "Frame stride", 1, 100, int(params.get("stride", 1)),
            key="wiz_stride2")

elif method in ("ABFE", "RBFE"):
    # OpenFE is the only ABFE/RBFE engine in the GUI, so every control here is
    # an OpenFE / OpenMM setting. Labels/help describe OpenFE behaviour.
    with st.expander("Sampling settings (OpenFE / OpenMM)", expanded=True):
        s1, s2, s3 = st.columns(3)
        params["complex_ns"] = s1.number_input(
            "Production per window (ns)", 0.1, 50.0, 5.0, 0.1, key="wiz_cns",
            help="OpenFE per-lambda-window production MD length "
                 "(HREX/repex). Complex and solvent legs use this value. "
                 "5 ns/window is a sane default; longer improves convergence.")
        params["abfe_equil_ns"] = s2.number_input(
            "Equilibration per window (ns)", 0.0, 10.0, 1.0, 0.1,
            key="wiz_abfe_eqns",
            help="OpenFE per-window equilibration before production. "
                 "Default 1 ns.")
        params["repeats"] = s3.number_input(
            "Repeats", 1, 5, int(_SITE.scheduler.max_concurrent > 0) or 3,
            key="wiz_rep",
            help="Independent repeats per transformation (for uncertainty).")
        # solvent leg mirrors the complex per-window length for the estimate.
        params["solvent_ns"] = params["complex_ns"]
        o1, o2 = st.columns(2)
        params["charge_method"] = o1.selectbox(
            "Ligand charge method", ["am1bcc", "nagl", "am1bccelf10"], 0,
            key="wiz_charge",
            help="Partial-charge scheme for the small molecule (OpenFE).")
        params["forcefield"] = o2.text_input(
            "Small-molecule force field", _SITE.openfe.small_molecule_ff,
            key="wiz_ff", help="OpenFF small-molecule force field, "
                               "e.g. openff-2.1.0.")
        if method == "RBFE":
            m1, m2 = st.columns(2)
            params["mapper"] = m1.selectbox(
                "Atom mapper", ["KartografAtomMapper", "LomapAtomMapper"],
                0, key="wiz_mapper",
                help="Maps common atoms across each ligand pair (RBFE).")
            params["network"] = m2.selectbox(
                "Network topology", ["generate_minimal_spanning_network",
                                     "generate_radial_network",
                                     "generate_maximal_network"], 0,
                key="wiz_net",
                help="How transformation edges are chosen between ligands.")

with st.expander("HPC / scheduler", expanded=False):
    st.caption("Defaults come from the **Settings** page (site config). "
               "Override here for this run only.")
    h1, h2 = st.columns(2)
    params["queue"] = h1.text_input("Queue / partition",
                                    _SITE.scheduler.gpu_queue, key="wiz_queue")
    params["walltime"] = h2.text_input("Walltime", _SITE.scheduler.walltime,
                                       key="wiz_wall")
    params["project"] = st.text_input(
        "Project / account", _SITE.scheduler.project, key="wiz_project",
        help="Accounting project/account passed to the scheduler "
             "(bsub -P / sbatch -A). Applies to all queue submissions.")
    # No client-side throttle here; the queue dispatches as GPUs free up.
    params["max_concurrent"] = 0
    st.caption("All jobs are submitted at once; the scheduler dispatches "
               "them as GPUs free up (see the per-page throttle for large runs).")
    if engine.startswith("OpenMM"):
        _ofe_py = _SITE.openfe.python_bin or (
            f"~/miniforge3/envs/{_SITE.openfe.conda_env}/bin/python")
        _ofe_bin = f"~/miniforge3/envs/{_SITE.openfe.conda_env}/bin/openfe"
        params["openfe_bin"] = st.text_input(
            "OpenFE `openfe` binary", _ofe_bin, key="wiz_ofebin")
        params["ofe_python"] = st.text_input(
            "OpenFE env python (for ABFE planner)", _ofe_py, key="wiz_ofepy")

# ============================ 4. Review & Launch ============================
st.subheader("4 - Review & launch")

ready = bool(protein and ligands and n_lig > 0 and ok)
if not ready:
    _missing = []
    if not protein: _missing.append("a protein")
    if not ligands or n_lig == 0: _missing.append("a ligand file")
    if not ok: _missing.append("a supported method/engine combo")
    if _missing:
        st.caption("To enable **Launch**, add: " + ", ".join(_missing) + ".")
if method == "RBFE" and eff_scope == "single":
    ready = False

# cost hint (very rough)
if method in ("ABFE", "RBFE"):
    units = (n_lig if method == "ABFE" else max(n_lig - 1, 0))
    jobs = units * int(params.get("repeats", 1))
    # final43: include MM-GBSA MD time (prod_ns) in the compute estimate;
    # complex_ns/solvent_ns only exist for ABFE/RBFE.
    if method == "MM-GBSA":
        total_ns = jobs * (params.get("prod_ns", 0) + params.get("equil_ns", 0))
    else:
        total_ns = jobs * (params.get("complex_ns", 0) + params.get("solvent_ns", 0))
    st.caption(f"Estimated: **{jobs} GPU job(s)**, ~{total_ns:.0f} ns total "
               f"sampling. All jobs submitted at once (LSF queue handles dispatch).")
elif method == "MM-GBSA":
    st.caption(f"Estimated: **{n_lig} MM-GBSA job(s)** ({engine}).")

rc1, rc2 = st.columns([1, 1])
validate = rc1.button("Validate", key="wiz_validate")
launch = rc2.button("Launch", type="primary", disabled=not ready, key="wiz_launch")

if validate:
    problems = []
    if not protein: problems.append("No protein selected.")
    if not ligands: problems.append("No ligand file selected.")
    if ligands and n_lig == 0: problems.append("Ligand file parsed to 0 records.")
    if not ok: problems.append(f"{method}/{engine} not supported.")
    if method == "RBFE" and eff_scope == "single":
        problems.append("RBFE requires >=2 ligands.")
    if ok and method == "MM-GBSA" and engine.startswith("OpenMM"):
        _vpf = preflight_openmm_mmgbsa(
            probe_python=str(Path(params.get("ofe_python", "python")).expanduser()),
            platform=str(params.get("platform", "CUDA")),
        )
        with st.expander("OpenMM MM-GBSA environment preflight",
                         expanded=not _vpf["ok"]):
            for k, v in _vpf["info"].items():
                st.write(f"**{k}**: `{v}`")
            for w in _vpf["warnings"]:
                st.warning(w)
            for er in _vpf["errors"]:
                st.error(er)
        if not _vpf["ok"]:
            problems.append("OpenMM MM-GBSA environment preflight failed "
                            "(see expander above).")
    if problems:
        for p in problems: st.error(p)
    else:
        st.success("All checks passed. Ready to launch.")


def _protonation_flags(P: dict) -> list[str]:
    """CLI protonation flags (config loader drops these fields, so they must
    go on the command line). --no-protonation disables auto-detect;
    --protonate CHAIN:RESNUM:NAME (repeatable) sets manual overrides."""
    flags: list[str] = []
    if not P.get("auto_protonation", True):
        flags.append("--no-protonation")
    raw = P.get("protonation_overrides", "") or ""
    for tok in raw.replace(",", "\n").splitlines():
        tok = tok.strip()
        if tok:
            flags += ["--protonate", tok]
    return flags


def _build_commands(exp: Experiment, wd: Path) -> list[list[str]]:
    """Translate the Experiment into one or more shell commands, reusing the
    existing engine entry points. Returns a list of argv lists."""
    # wd is the already-timestamped, method-specific run directory supplied by
    # the caller (see the launch block). All --workdir/--out args derive from it.
    wd.mkdir(parents=True, exist_ok=True)
    P = exp.params
    cmds: list[list[str]] = []

    if exp.method == "RBFE" and exp.engine.startswith("OpenMM"):
        # 1) plan the network (login-node). Run/quickrun happens on Monitor page.
        net = wd / "network_setup"
        yaml_path = write_network_yaml(
            OpenFESettings(
                sim_time_ns=P.get("complex_ns", 10.0),
                equil_time_ns=1.0, n_replicates=int(P.get("repeats", 1)),
                sampler="repex", platform="CUDA",
                forcefield=P.get("forcefield", "openff-2.2.0"),
                mapper=P.get("mapper", "KartografAtomMapper"),
                network=P.get("network", "generate_minimal_spanning_network"),
                charge_method=P.get("charge_method", "am1bcc"),
                project=P.get("project", "your-project"),
                queue=P.get("queue", "gpu"),
                walltime=P.get("walltime", "48:00"), conda_env="openfe_env",
                max_concurrent=int(P.get("max_concurrent", 0)),
            ),
            wd / "openfe_network.yaml")
        cmds.append(plan_rbfe_cmd(
            Path(exp.protein), Path(exp.ligands), net, yaml_path,
            openfe_bin=str(Path(P.get("openfe_bin", "openfe")).expanduser())))

    elif exp.method == "RBFE" and exp.engine == "Amber":
        # Removed from the GUI (v2.6.0). Amber RBFE is not an exposed,
        # supported workflow; the engine code remains in amber_md for CLI use.
        raise RuntimeError(
            "Amber RBFE is not available from the GUI. Use RBFE / OpenMM / OpenFE.")

    elif exp.method == "ABFE" and exp.engine.startswith("OpenMM"):
        # OpenFE absolute-binding planner (runs in the OpenFE env python).
        py = str(Path(P.get("ofe_python", "python")).expanduser())
        cmds.append([
            py, "-m", "amber_md.abfe_openfe_plan",
            "--protein", str(Path(exp.protein)),
            "--ligands", str(Path(exp.ligands)),
            "--out", str(wd / "abfe_setup"),
            "--complex-ns", str(P.get("complex_ns", 10.0)),
            "--solvent-ns", str(P.get("solvent_ns", 5.0)),
            "--repeats", str(int(P.get("repeats", 1))),
            "--charge-method", P.get("charge_method", "am1bcc"),
            "--forcefield", P.get("forcefield", "openff-2.2.0"),
        ])

    elif exp.method == "MM-GBSA" and exp.engine.startswith("OpenMM"):
        # OpenMM MD + AmberTools MMPBSA.py (Interpretation A).
        # Multi-ligand fan-out (Option 2): split the ligand file into per-ligand
        # single-record SDFs and emit ONE job per ligand into
        #   <work_dir>/lig_<name>/   (matches batch_aggregate's lig_*/mmgbsa/...)
        # so the existing ranking (amber_md.batch_aggregate) works unchanged.
        # Single-ligand inputs naturally produce a single job.
        from amber_md.mmgbsa_openmm import split_ligand_file
        py = str(Path(P.get("ofe_python", "python")).expanduser())
        try:
            recs = split_ligand_file(Path(exp.ligands), wd / "_ligands")
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(f"Could not split ligand file: {e}")
        common = [
            # final43: read the explicit MM-GBSA time fields first; fall back
            # to complex_ns (legacy) then the 5 ns / 1 ns hard default.
            "--prod-ns", str(P.get("prod_ns", P.get("complex_ns", 10.0))),  # final51: unified 10 ns
            "--equil-ns", str(P.get("equil_ns", 1.0)),
            "--igb", str(P.get("igb", 8)),
            "--salt", str(P.get("salt_conc", 0.15)),
            "--platform", str(P.get("platform", "CUDA")),
            "--buffer", str(P.get("buffer_A", 10.0)),
            "--charge-method", str(P.get("mm_charge_method", "bcc")),
            "--ligand-charge", str(int(P.get("ligand_charge", 0))),
        ]
        if P.get("decomposition"):
            common.append("--decomp")
        # Execution target -> LSF GPU queue submission (per-ligand job is
        # submitted by its own CLI invocation, which writes a #BSUB script and
        # calls bsub, then exits). Default 'Local host' keeps the old behaviour.
        if P.get("exec_target", "Local host").startswith("GPU queue"):
            common += ["--submit", "gpu",
                       "--queue", str(P.get("gpu_queue", "gpu")),
                       "--walltime", str(P.get("gpu_walltime", "24:00")),
                       "--conda-env", str(Path(P.get("ofe_python", "python"))
                                          .expanduser().parent.parent)]
            _proj = str(P.get("project", "your-project")).strip()
            if _proj:
                common += ["--project", _proj]
        for r in recs:
            lig_dir = wd / f"lig_{r['name']}"
            c = [
                py, "-m", "amber_md.mmgbsa_openmm",
                "--protein", str(Path(exp.protein)),
                "--ligand", str(r["sdf"]),
                "--workdir", str(lig_dir),
            ] + common
            cmds.append(c)

    elif exp.method == "MM-GBSA" and exp.engine == "Amber":
        # run_amber.py MM-GBSA, FANNED OUT one job per ligand (final59).
        # Previously the whole multi-record ligand file went to a single
        # invocation -> one job (first record only). MM-GBSA must submit one
        # job per molecule. Reuse the OpenMM path's splitter; emit one command
        # per ligand into <wd>/lig_<name>/ (matches batch_aggregate layout).
        from amber_md.mmgbsa_openmm import split_ligand_file
        try:
            recs = split_ligand_file(Path(exp.ligands), wd / "_ligands")
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(f"Could not split ligand file: {e}")
        for r in recs:
            lig_dir = wd / f"lig_{r['name']}"
            cfg_path = amber_config.write_config(
                P, work_dir=str(lig_dir), method="MM-GBSA",
                protein=str(Path(exp.protein)), ligand=str(r["sdf"]))
            c = [
                "python", str(run_amber_py()),
                "--config", str(cfg_path),
                "--protein-file", str(Path(exp.protein)),
                "--ligand-file", str(r["sdf"]),
                "--workdir", str(lig_dir),
                "--salt", str(P.get("salt_conc", 0.15)),
                "--prod-ns", str(P.get("prod_ns", 10.0)),
                "--equil-ns", str(P.get("equil_ns", 1.0)),
            ]
            _rn_mm = str(P.get("ligand_resname", "")).strip()
            if _rn_mm:
                c += ["--lig-resname", _rn_mm]
            if P.get("decomposition"):
                c.append("--decomp")
            c += _protonation_flags(P)
            if P.get("queue"):
                c += ["--queue", str(P["queue"])]
            if P.get("walltime"):
                c += ["--walltime", str(P["walltime"])]
            cmds.append(c)

    elif exp.method == "ABFE" and exp.engine == "Amber":
        # Removed from the GUI (v2.6.0). Amber ABFE is not an exposed,
        # supported workflow; the engine code remains in amber_md for CLI use.
        raise RuntimeError(
            "Amber ABFE is not available from the GUI. Use ABFE / OpenMM / OpenFE.")

    return cmds


if launch:
    # Surface the execution-target choice into params so _build_commands can
    # translate it into --submit/--queue/--walltime for OpenMM MM-GBSA.
    params = dict(params)
    params["exec_target"] = exec_target
    params["gpu_queue"] = gpu_queue
    params["gpu_walltime"] = gpu_walltime
    params["throttle_on"] = bool(throttle_on)
    params["throttle_n"] = int(throttle_n)
    # When GPU queue is chosen, route the existing queue-aware paths (Amber
    # MM-GBSA/ABFE, OpenFE RBFE) through the selected queue/walltime WITHOUT
    # changing their behaviour when the user leaves Local host selected.
    if exec_target.startswith("GPU queue"):
        params.setdefault("queue", gpu_queue)
        params.setdefault("walltime", gpu_walltime)
    exp = Experiment(
        method=method, engine=engine, scope=scope_choice,
        protein=protein, ligands=ligands, work_dir=work_dir,
        n_ligands=n_lig, effective_scope=eff_scope, params=params,
    )
    # Per-run isolation: each method gets its own <method>_<timestamp> directory
    # under the base work_dir, so MM-GBSA / RBFE / ABFE runs never share a folder
    # (even multiple launches on the same day).
    _method_tag = {
        "MM-GBSA": "mmgbsa", "RBFE": "rbfe", "ABFE": "abfe",
    }.get(exp.method, exp.method.lower().replace("/", "_").replace(" ", "_"))
    run_wd = (Path(work_dir).expanduser()
              / f"{_method_tag}_{_dt.now().strftime('%Y%m%d_%H%M%S')}")
    try:
        cmds = _build_commands(exp, run_wd)
    except Exception as e:                       # noqa: BLE001
        st.error(f"Failed to build launch command: {e}")
        cmds = []

    # ---- Preflight: OpenFE RBFE AM1BCC charge generation needs a working
    # AmberTools (antechamber) inside the OpenFE env. Catch the common
    # "AmberTools version None" / broken system-Amber failure BEFORE spawning.
    if cmds and exp.method == "RBFE" and exp.engine.startswith("OpenMM"):
        _pf = preflight_openfe_charges(
            openfe_bin=str(Path(params.get("openfe_bin", "openfe")).expanduser()),
            charge_method=params.get("charge_method", "am1bcc"),
        )
        with st.expander("OpenFE charge-generation preflight", expanded=not _pf["ok"]):
            for k, v in _pf["info"].items():
                st.write(f"**{k}**: `{v}`")
            for w in _pf["warnings"]:
                st.warning(w)
            for er in _pf["errors"]:
                st.error(er)
        if not _pf["ok"]:
            st.error(
                "Preflight failed: the OpenFE environment cannot generate "
                "AM1BCC partial charges, so `openfe plan-rbfe-network` would "
                "crash in antechamber. Fix the env (see messages above) -- "
                "typically:  conda install -n <openfe_env> -c conda-forge "
                "ambertools  and `unset AMBERHOME` -- or pick an RDKit charge "
                "method (gasteiger) to bypass AmberTools. Launch aborted.")
            cmds = []

    # ---- Preflight: OpenMM MM-GBSA needs OpenMM (+ CUDA platform) AND the
    # AmberTools scoring binaries in the same env. Catch missing OpenMM / CUDA /
    # AmberTools BEFORE spawning the (long) MD job.
    if cmds and exp.method == "MM-GBSA" and exp.engine.startswith("OpenMM"):
        _pf = preflight_openmm_mmgbsa(
            probe_python=str(Path(params.get("ofe_python", "python")).expanduser()),
            platform=str(params.get("platform", "CUDA")),
        )
        with st.expander("OpenMM MM-GBSA preflight", expanded=not _pf["ok"]):
            for k, v in _pf["info"].items():
                st.write(f"**{k}**: `{v}`")
            for w in _pf["warnings"]:
                st.warning(w)
            for er in _pf["errors"]:
                st.error(er)
        if not _pf["ok"]:
            st.error(
                "Preflight failed: the env cannot run OpenMM MM-GBSA (missing "
                "OpenMM, the requested platform, or AmberTools). Fix the env "
                "(see messages above) -- typically: conda install -n <env> "
                "-c conda-forge openmm ambertools mdtraj -- then retry. "
                "Launch aborted.")
            cmds = []

    if not cmds:
        st.error("No command was generated for this combination yet.")
    else:
        wd = run_wd                       # the timestamped, method-specific dir
        wd.mkdir(parents=True, exist_ok=True)
        st.info(f"This run's working directory: `{wd}`")
        launched = []
        # For OpenFE RBFE, spawn with a sanitized env so the detached process
        # uses the conda env's antechamber, not a system Amber inherited from
        # the launching shell (the "AmberTools version None" failure mode).
        spawn_env = None
        if exp.method == "RBFE" and exp.engine.startswith("OpenMM"):
            spawn_env = sanitized_openfe_env(
                str(Path(params.get("openfe_bin", "openfe")).expanduser()))
        elif exp.method == "MM-GBSA" and exp.engine.startswith("OpenMM"):
            # OpenMM MM-GBSA runs in the OpenFE env (OpenMM + AmberTools);
            # sanitize so a system Amber inherited from the shell does not shadow
            # the env's tleap/antechamber/MMPBSA.py.
            spawn_env = sanitized_openfe_env(
                str(Path(params.get("ofe_python", "python")).expanduser()))
        def _cmd_label(cmd, i):
            # Use the --workdir/--work-dir basename (e.g. lig_<name>) when present
            # so each fan-out job has a readable, unique log name. final59: also
            # match the fep_driver spelling "--work-dir" (ABFE/RBFE) so per-ligand
            # ABFE jobs no longer collide on a single launch_abfe_0.log.
            for flag in ("--workdir", "--work-dir"):
                try:
                    wi = cmd.index(flag)
                    base = Path(cmd[wi + 1]).name
                    if base:
                        return base
                except (ValueError, IndexError):
                    pass
            return f"{exp.method.lower()}_{i}"
        _is_queue = (
            params.get("exec_target", "Local host").startswith("GPU queue")
            and (
                (exp.method == "MM-GBSA" and exp.engine.startswith("OpenMM"))
                # final56: ABFE/Amber now submits to LSF too (driver --submit).
                or (exp.method == "ABFE" and exp.engine == "Amber")
                # final58: RBFE/Amber per-edge legs submit via fep_driver too.
                or (exp.method == "RBFE" and exp.engine == "Amber")
            ))
        for i, cmd in enumerate(cmds):
            log = wd / f"launch_{_cmd_label(cmd, i)}.log"
            st.code(" ".join(shlex.quote(c) for c in cmd), language="bash")
            if exp.method in ("ABFE", "RBFE") and exp.engine == "Amber":
                st.caption("✓ Pre-submit mdin validation runs automatically "
                           "in the driver (catches malformed heat/dens/eq cards "
                           "before any GPU job). Check the launch log for its result.")
            try:
                pid = spawn_detached(cmd, log, cwd=str(wd), env=spawn_env)
                launched.append((pid, str(log)))
                if _is_queue:
                    st.success(f"Submitting to LSF GPU queue "
                               f"'{params.get('gpu_queue','gpu')}' "
                               f"(bsub launcher pid {pid}). Job-id in: {log}")
                else:
                    st.success(f"Launched (pid {pid}). Log: {log}")
            except Exception as e:               # noqa: BLE001
                st.error(f"Launch failed: {e}")
        if _is_queue:
            st.info("Jobs were submitted to LSF. Track them with `bjobs` or the "
                    "Job Monitor page. Each per-ligand job writes its own "
                    "`lsf_mmgbsa.<jobid>.out/.err` in its lig_* directory.")
        if exp.method == "MM-GBSA" and exp.engine.startswith("OpenMM") and len(cmds) > 1:
            st.info(
                f"Fanned out **{len(cmds)} per-ligand MM-GBSA jobs** into "
                f"`{wd}/lig_*/`. When they finish, rank them with:  "
                f"`python -m amber_md.batch_aggregate {wd}`  "
                "(opens INDEX.html with the sorted binding-energy table). "
                "Each per-ligand result lands at "
                "`lig_<name>/mmgbsa/FINAL_RESULTS_MMPBSA.dat`.")
        # persist for the Monitor / Results pages
        st.session_state["experiment"] = asdict(exp)
        _ll = {
            "pids": [p for p, _ in launched],
            "logs": [l for _, l in launched],
            "work_dir": str(wd),
            "method": exp.method, "engine": exp.engine,
            "scope": exp.effective_scope,
            "ts": time.time(),
        }
        st.session_state["last_launch"] = _ll
        # durable registry consumed by the Job Monitor page (Phase 2)
        st.session_state.setdefault("launches", []).append(_ll)
        if launched:
            st.info("Open **Job Monitor** to track progress.")

# footer: quick glance at active LSF jobs (reuses existing helper)
st.divider()
with st.expander("Active LSF jobs (quick glance)", expanded=False):
    jobs = get_lsf_jobs()
    if jobs:
        st.dataframe(jobs, hide_index=True, width="stretch")
    else:
        st.caption("No active LSF jobs.")
