"""
5_Results_Single.py  --  Results: Single Molecule (v2.5.0, Phase 3).

Detailed per-system analysis for ONE ligand/run:
  * Headline DG_bind (with honest uncertainty)
  * Energetic breakdown (per-leg: complex / solvent / SSC)
  * Convergence diagnostics (MBAR overlap; fwd/rev if present)
  * MM-GBSA per-residue decomposition (when present)

Auto-targets the ligand recorded by the Setup & Launch wizard
(st.session_state["experiment"]); a directory picker overrides it.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="Results: Single", layout="wide",
                   page_icon="microscope")

from amber_md.gui.common import dir_picker
from amber_md.gui import results_lib as rl
from amber_md.gui import openfe_report as ofr
from amber_md import mmpbsa_report as mpr  # final50: unified report button


st.title("Results - Single Molecule")
st.caption("Detailed energetics and convergence for one run. "
           "Pick a run directory, or it follows your last launch.")

# ---- resolve target dir: session experiment -> picker ----
exp = st.session_state.get("experiment", {})
default_wd = exp.get("work_dir", "~/Run_dir")
with st.sidebar:
    st.header("Run directory")
    wd = dir_picker("Run / ligand directory", "rs_dir", default_dir=default_wd)
wd = Path(wd).expanduser() if wd else None

if not wd or not wd.exists():
    st.info("Select a run directory in the sidebar to view its results.")
    st.stop()

st.subheader(f"Run: `{wd.name}`")
if exp:
    st.caption(f"Context from last launch: **{exp.get('method')}** / "
               f"**{exp.get('engine')}** ({exp.get('effective_scope')})")

# ---- locate result source ----
fep_dir = wd / "fep"
mm_dat = wd / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"

# Replicate-aware OpenFE discovery: handles both <ligand>_rep*/results.json and
# flat *_result.json. If the picked dir holds several replicates of one ligand
# (or several ligands), let the user choose which run to view in detail.
_groups = ofr.discover_ligands(wd)
openfe_candidates = []          # ordered list of Paths, selected one is [0]
_active_reps = []               # ReplicateResult list for the chosen ligand
if _groups:
    # pick ligand (default: the one matching wd.name, else first)
    _lig_keys = sorted(_groups.keys())
    _default_key = next((k for k in _lig_keys if k in wd.name or wd.name in k),
                        _lig_keys[0])
    if len(_lig_keys) > 1:
        _sel_lig = st.sidebar.selectbox("Ligand", _lig_keys,
                                        index=_lig_keys.index(_default_key),
                                        key="rs_lig")
    else:
        _sel_lig = _default_key
    _active_reps = _groups[_sel_lig]
    # choose replicate
    _rep_labels = [r.label for r in _active_reps]
    if len(_rep_labels) > 1:
        _sel_rep = st.sidebar.selectbox(
            f"Replicate ({len(_rep_labels)} found)", _rep_labels, index=0,
            key="rs_rep")
        # mean +/- SEM banner across all replicates of this ligand
        _st = ofr.replicate_stats(_active_reps)
        if _st["mean"] is not None:
            if _st["sem"] is not None:
                st.info(f"**{_sel_lig}**: {len(_rep_labels)} replicates -> "
                        f"mean DG_bind **{_st['mean']:.2f} +/- {_st['sem']:.2f}** "
                        f"kcal/mol (SD {_st['sd']:.2f}). Showing replicate "
                        f"`{_sel_rep}` below; use the Compare page or Build "
                        f"report for the full multi-replicate view.")
            else:
                st.info(f"**{_sel_lig}**: showing replicate `{_sel_rep}`.")
    else:
        _sel_rep = _rep_labels[0] if _rep_labels else None
    # resolve the chosen replicate's result path to the front of the list
    _chosen = next((r for r in _active_reps if r.label == _sel_rep),
                   _active_reps[0] if _active_reps else None)
    if _chosen is not None:
        from pathlib import Path as _P
        openfe_candidates = [_P(_chosen.path)] + [
            _P(r.path) for r in _active_reps if r.path != _chosen.path]

tab_summary, tab_breakdown, tab_conv, tab_decomp = st.tabs(
    ["Headline", "Energy breakdown", "Convergence", "MM-GBSA decomp"])

# =========================== HEADLINE ===========================
with tab_summary:
    shown = False
    # OpenFE
    if openfe_candidates:
        rj = openfe_candidates[0]
        r = rl.parse_openfe_result(rj)
        if r:
            shown = True
            est = r["estimate_kcal_mol"]
            unc = r["uncertainty_kcal_mol"]
            # honest error: propagate per-leg if reported uncertainty is 0
            import math
            legs = r["legs"]
            prop = None
            if legs:
                errs = [v.get("unit_estimate_error") or 0.0 for v in legs.values()]
                prop = math.sqrt(sum(e * e for e in errs)) if errs else None
            disp_err = unc if (unc and unc > 0) else prop
            c1, c2 = st.columns([2, 1])
            c1.metric("ΔG_bind (kcal/mol)",
                      f"{est:.2f}" if est is not None else "—",
                      f"± {disp_err:.2f}" if disp_err else None)
            if (not unc) and prop:
                c2.warning("Reported uncertainty is 0 (single repeat). "
                           f"Showing propagated per-leg error ±{prop:.2f}.")
            st.caption(f"Source: `{rj.name}`")
    # Amber FEP
    if not shown and fep_dir.exists():
        val, kind, legs = rl.fep_headline(fep_dir)
        if val is not None:
            shown = True
            st.metric(f"{kind} (kcal/mol)", f"{val:.2f}")
    # MM-GBSA  (final64: show the detected engine via the locked detector)
    if not shown and mm_dat.exists():
        status, dg, full = rl.mmgbsa_status(wd)
        if dg is not None:
            shown = True
            engine = rl.mmgbsa_engine(wd)
            c1, c2 = st.columns([2, 1])
            c1.metric("MM-GBSA ΔG_bind (kcal/mol)", f"{dg:.2f}",
                      f"± {full.get('sem'):.2f}" if full and full.get("sem")
                      else None)
            c2.metric("Engine", engine)
            st.caption(f"Method: **MM-GBSA** / Engine: **{engine}** "
                       "(engine.json marker, else OpenMM MD artifacts, "
                       "else Amber).")
    if not shown:
        st.info("No parseable result found in this directory yet.")

# =========================== BREAKDOWN ===========================
with tab_breakdown:
    if openfe_candidates:
        r = rl.parse_openfe_result(openfe_candidates[0])
        if r and r["legs"]:
            rows = []
            for leg, d in r["legs"].items():
                rows.append({"leg": leg,
                             "ΔG (kcal/mol)": d.get("unit_estimate"),
                             "± err": d.get("unit_estimate_error"),
                             "SSC": d.get("standard_state_correction"),
                             "prod iters": d.get("production_iterations")})
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    elif fep_dir.exists():
        legs = rl.fep_run_results(fep_dir)
        if legs:
            rows = [{"leg": k, "ΔG (kcal/mol)": v.get("dG_kcal_mol"),
                     "± err": v.get("err_kcal_mol"),
                     "estimator": v.get("estimator")} for k, v in legs.items()]
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            for k, v in legs.items():
                if "estimators_df" in v:
                    st.markdown(f"**{k} - estimator comparison**")
                    st.dataframe(v["estimators_df"], hide_index=True,
                                 width="stretch")
    else:
        st.info("No per-leg breakdown available.")

# =========================== CONVERGENCE ===========================
with tab_conv:
    if openfe_candidates:
        r = rl.parse_openfe_result(openfe_candidates[0])
        if r and r["legs"]:
            rows = []
            for leg, d in r["legs"].items():
                ov = d.get("mbar_overlap_scalar")
                rows.append({"leg": leg, "MBAR overlap (scalar)": ov,
                             "equil iters": d.get("equilibration_iterations"),
                             "prod iters": d.get("production_iterations"),
                             "health": ("OK" if (ov or 0) > 0.005 else "check")})
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            st.caption("Scalar overlap is a global metric; adjacent-window "
                       "overlap (from the full matrix) is the stricter check.")
        else:
            st.info("No convergence diagnostics in this result.")
    else:
        st.info("Convergence diagnostics are available for OpenFE results.")

# =========================== MM-GBSA DECOMP ===========================
with tab_decomp:
    decomp = wd / "mmgbsa" / "FINAL_DECOMP_MMPBSA.dat"
    if decomp.exists():
        st.caption(f"`{decomp}`")
        txt = decomp.read_text(errors="ignore")
        st.code(txt[:4000] + ("\n... (truncated)" if len(txt) > 4000 else ""))
    else:
        st.info("No FINAL_DECOMP_MMPBSA.dat (run MM-GBSA with decomposition "
                "enabled to populate this).")


# ---------------------------------------------------------------------------
# final50: ONE "Generate HTML report" button that auto-detects the result type
# (OpenFE campaign  vs  Amber/OpenMM MM-GBSA) and calls the right generator.
# ---------------------------------------------------------------------------
def _find_mmgbsa_dats(root):
    """Return MM-GBSA .dat files under `root`, tolerating MMPBSA/MMGBSA naming
    and single (mmgbsa/) vs batch (lig_*/mmgbsa/) and workdir-root layouts."""
    names = ("FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMGBSA.dat")
    hits = []
    for sub in ("mmgbsa", "."):
        for nm in names:
            p = root / sub / nm
            if p.exists() and p not in hits:
                hits.append(p)
    for nm in names:
        hits += [p for p in sorted(root.glob(f"lig_*/mmgbsa/{nm}")) if p not in hits]
    return hits


def _detect_report_kind(root):
    """('openfe', n_ligands) | ('mmgbsa', [dat,...]) | (None, None)."""
    try:
        groups = ofr.discover_ligands(root)
    except Exception:
        groups = {}
    if groups:
        return "openfe", len(groups)
    dats = _find_mmgbsa_dats(root)
    if dats:
        return "mmgbsa", dats
    return None, None

st.divider()
st.subheader("HTML report")

_kind, _info = _detect_report_kind(wd)
if _kind == "openfe":
    st.caption(f"Detected **OpenFE** results ({_info} ligand group(s)).")
elif _kind == "mmgbsa":
    st.caption(f"Detected **MM-GBSA** results ({len(_info)} .dat file(s)).")
    # surface any reports that already exist
    _existing = [d.with_name("FINAL_RESULTS.report.html") for d in _info
                 if d.with_name("FINAL_RESULTS.report.html").exists()]
    if _existing:
        st.info(f"{len(_existing)} report(s) already generated; click to "
                "rebuild, or download below.")
        for h in _existing:
            try:
                st.download_button(f"Download {h.parent.name}/{h.name}",
                                   h.read_text(), file_name=h.name,
                                   mime="text/html", key=f"rs_dl_exist_{h}")
            except Exception:
                pass
else:
    st.caption("No OpenFE or MM-GBSA results detected in this directory yet.")

if st.button("Generate HTML report", key="rs_build_report",
             disabled=(_kind is None)):
    try:
        if _kind == "openfe":
            with st.spinner("Discovering replicates and rendering OpenFE HTML..."):
                res = ofr.build_campaign_report(wd)
            rows = res["rows"]
            if not rows:
                st.warning("No *_result.json / results.json found under this directory.")
            else:
                st.success(f"OpenFE report written for {len(rows)} ligand(s). "
                           f"Index: {res['index']}")
                try:
                    idx_html = Path(res["index"]).read_text()
                    st.download_button("Download campaign index.html", idx_html,
                                       file_name="index.html", mime="text/html",
                                       key="rs_report_dl")
                except Exception:
                    pass
                st.info(f"Per-ligand HTML files are in: {Path(res['index']).parent}")
        elif _kind == "mmgbsa":
            built = []
            with st.spinner(f"Rendering MM-GBSA HTML for {len(_info)} result(s)..."):
                for dat in _info:
                    try:
                        built.append(mpr.generate_report(dat))
                    except Exception as e:  # noqa: BLE001
                        st.warning(f"Skipped {dat}: {e}")
            if built:
                st.success(f"MM-GBSA report(s) written: {len(built)} file(s).")
                for h in built:
                    try:
                        st.download_button(
                            f"Download {Path(h).parent.name}/{Path(h).name}",
                            Path(h).read_text(), file_name=Path(h).name,
                            mime="text/html", key=f"rs_dl_{h}")
                    except Exception:
                        pass
            else:
                st.error("No MM-GBSA report could be generated (see warnings).")
    except Exception as e:  # noqa: BLE001
        st.error(f"Report generation failed: {e}")
