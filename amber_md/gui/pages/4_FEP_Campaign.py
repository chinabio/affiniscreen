"""
4_FEP_Campaign.py  --  OpenFE FEP Campaign driver (RBFE / ABFE).

Drives the OpenFE / OpenMM free-energy campaign end to end:
  run -> gather -> solve -> cycle-closure.

The one bsub-per-(edge x replicate) submission is handled by
``amber_md.gui.openfe_campaign``; solve + cycle-closure are engine-agnostic
(``amber_md.gui.fep_common``).

Note: the legacy Amber TI engine selector was removed with the rest of the
Amber ABFE/RBFE GUI surface. The Amber engine code is still shipped in the
package for programmatic use.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Part of AffiniScreen.

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="FEP Campaign", layout="wide", page_icon="rocket")

from amber_md.gui.common import dir_picker, spawn_detached
from amber_md.gui.openfe_common import OpenFESettings, MAX_GPU
from amber_md.gui.fep_common import solve_network, cycle_closure_residuals
from amber_md.gui import openfe_campaign as oc_openfe
from amber_md import site_config as _sc
_SITE = _sc.get()


st.title("OpenFE FEP Campaign - Run, Gather, Solve (RBFE / ABFE)")
st.caption("OpenFE / OpenMM free-energy campaign: "
           "run -> gather -> solve -> cycle-closure.")

exp = st.session_state.get("experiment", {})
default_wd = exp.get("work_dir", "~/Run_dir/run_v250")

# OpenFE / OpenMM is the only FEP-campaign engine exposed in the GUI.
engine = "OpenFE / OpenMM"
oc = oc_openfe
with st.sidebar:
    st.header("Engine")
    st.info("OpenFE / OpenMM")

    st.divider()
    st.header("Campaign directory")
    wd = dir_picker("Work directory (contains network_setup/ or abfe_setup/)",
                    "oc_wd", default_dir=default_wd)
    st.divider()
    st.header("Run settings")
    n_rep = st.number_input("Replicates per edge", 1, 5,
                            int(exp.get("params", {}).get("repeats", 3)),
                            key="oc_nrep")
    queue = st.text_input("Queue / partition", _SITE.scheduler.gpu_queue,
                          key="oc_queue")
    project = st.text_input("Project / account", _SITE.scheduler.project,
                            key="oc_proj")
    walltime = st.text_input("Walltime", _SITE.scheduler.walltime, key="oc_wall")
    conda_env = st.text_input("OpenFE conda env", _SITE.openfe.conda_env,
                              key="oc_env")

wd = Path(wd).expanduser() if wd else None
if not wd or not wd.exists():
    st.info("Select a campaign work directory (the one the wizard planned into).")
    st.stop()

layout = oc.CampaignLayout.at(wd)
settings = OpenFESettings(
    sim_time_ns=exp.get("params", {}).get("complex_ns", 10.0),
    equil_time_ns=1.0, n_replicates=int(n_rep), sampler="repex",
    platform="CUDA", forcefield=_SITE.openfe.small_molecule_ff,
    mapper="KartografAtomMapper",
    network="generate_minimal_spanning_network", charge_method="am1bcc",
    project=project, queue=queue, walltime=walltime, conda_env=conda_env,
    extra_modules=[], max_concurrent=0)

es = oc.edges(layout)
if not es:
    st.warning(f"No edges found for **{engine}** under `{layout.network_dir}`. "
               "Plan first via Setup & Launch (writes network_setup/ or "
               "abfe_setup/).")
    st.stop()

st.success(f"**{engine}** - found **{len(es)}** edge(s) in "
           f"`{layout.network_dir.name}`.")

tab_map, tab_run, tab_status, tab_solve, tab_script = st.tabs(
    ["Atom mapping", "Run edges", "Status", "Gather & solve", "Preview script"])

with tab_map:
    st.markdown("#### Atom-mapping inspection")
    arts = oc.find_mapping_artifacts(layout)
    if not arts:
        st.info("No mapping artifacts found yet. Plan the network first.")
    else:
        if "diagnostics" in arts:
            with st.expander("Network diagnostics", expanded=True):
                st.code(arts["diagnostics"].read_text(errors="ignore"))
        if "edges" in arts:
            rows = oc.read_edges_table(arts["edges"])
            st.markdown(f"**{len(rows)} edge(s)**")
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            st.markdown("##### Perturbation network")
            st.markdown("```mermaid\n" + oc.mapping_network_mermaid(rows) + "\n```")
            with st.expander("Per-edge MCS SMARTS & masks"):
                for e in rows:
                    st.markdown(f"**{e.get('lig_a','?')} ~ {e.get('lig_b','?')}** "
                                f"- core_heavy={e.get('core_heavy','')}, "
                                f"perturbed={e.get('perturbed_atoms','')}, "
                                f"score={e.get('score','')}")
                    if e.get("mcs_smarts"):
                        st.code(e["mcs_smarts"], language="text")
        else:
            st.warning("Found edges but no edges.csv with mapping diagnostics.")
        if "transformations_dir" in arts:
            with st.expander(f"Edges ({len(es)})"):
                st.dataframe(pd.DataFrame([{"edge": e.name} for e in es]),
                             hide_index=True, width="stretch")

with tab_run:
    n_jobs = len(es) * int(n_rep)
    st.markdown(f"**{len(es)} edges x {n_rep} replicate(s) = {n_jobs} GPU jobs.**")
    st.caption(f"Currently {oc.count_my_jobs(oc.JOB_PREFIX)} "
               f"`{oc.JOB_PREFIX}*` jobs in the queue.")
    _tc1, _tc2 = st.columns([1.3, 1])
    oc_throttle_on = _tc1.checkbox("Throttle submissions", value=False,
                                   key="oc_throttle_on")
    oc_throttle_n = _tc2.number_input("Max jobs in queue", 1, 200, 8, 1,
                                      key="oc_throttle_n",
                                      disabled=not oc_throttle_on)
    if st.button("Submit all (edge x repeat)", type="primary",
                 key="oc_submit"):
        progress = st.progress(0.0); status = st.empty()

        def _cb(i, total, sub, skip):
            progress.progress(i / total)
            status.write(f"Submitted {sub} / {total} (skipped {skip} done)")

        sub, skip, errs = oc.submit_campaign(
            layout, settings, int(n_rep), progress_cb=_cb,
            throttle_n=(int(oc_throttle_n) if oc_throttle_on else 0))
        if oc_throttle_on and len(errs) == 1 and errs[0].startswith("THROTTLED"):
            st.success(f"Handed {sub} job(s) to the login-node throttler.")
            st.info(errs[0]); errs = []
        for e in errs:
            st.error(f"Submission failed - {e}")
        if not (oc_throttle_on and sub):
            st.success(f"Done. Submitted {sub}, skipped {skip} complete.")
        st.session_state.setdefault("launches", []).append(
            {"method": "RBFE", "engine": "OpenMM / OpenFE",
             "work_dir": str(wd), "pids": [], "logs": [], "ts": time.time()})

with tab_status:
    auto = st.checkbox("Auto-refresh every 30 s", value=False, key="oc_auto")
    snap = oc.campaign_status(layout, int(n_rep))
    st.progress(snap["pct"], text=f"DONE {snap['done']} / {snap['total']} "
                                   f"({100*snap['pct']:.0f}%)")
    b = snap["buckets"]; c = st.columns(4)
    c[0].metric("DONE", b.get("DONE", 0))
    c[1].metric("analyzing", b.get("analyzing", 0))
    c[2].metric("running", b.get("running", 0))
    c[3].metric("not started", b.get("not started", 0))
    st.dataframe(pd.DataFrame(snap["rows"]), hide_index=True, width="stretch")
    if auto:
        time.sleep(30); st.rerun()

with tab_solve:
    st.markdown("#### 1) Gather per-edge results")
    cmd = oc.gather_command(layout)
    st.code(" ".join(str(x) for x in cmd), language="bash")
    if st.button("Run `openfe gather`", key="oc_gather"):
        log = layout.work_dir / "gather.log"
        spawn_detached(cmd, log, cwd=str(layout.work_dir))
        st.success(f"Gather launched -> `{layout.dg_tsv}` (log: {log}).")

    st.markdown("#### 2) Edge results table")
    edge_tsv = layout.dg_tsv
    edges_with_ddg = []
    if edge_tsv.exists():
        try:
            df = pd.read_csv(edge_tsv, sep="\t")
            st.dataframe(df, hide_index=True, width="stretch")
            cols = {c.lower(): c for c in df.columns}
            a_c = next((cols[k] for k in cols if k in
                        ("liga", "ligand_i", "ligand_a", "from")), None)
            b_c = next((cols[k] for k in cols if k in
                        ("ligb", "ligand_j", "ligand_b", "to")), None)
            dg_c = next((cols[k] for k in cols if "ddg" in k or
                         k in ("dg(kcal/mol)", "estimate")), None)
            er_c = next((cols[k] for k in cols if "unc" in k or "err" in k or
                         "stderr" in k), None)
            if a_c and b_c and dg_c:
                for _, r in df.iterrows():
                    try:
                        edges_with_ddg.append((
                            str(r[a_c]), str(r[b_c]), float(r[dg_c]),
                            float(r[er_c]) if er_c and pd.notna(r[er_c]) else 0.5))
                    except (ValueError, TypeError):
                        continue
        except Exception as e:  # noqa: BLE001
            st.warning(f"Could not parse {edge_tsv}: {e}")
    else:
        st.info("Run gather first to produce the edge TSV.")

    st.markdown("#### 3) Solve network -> per-ligand DeltaG")
    if edges_with_ddg:
        ref = st.text_input("Reference ligand (pinned to 0; blank = first)",
                            "", key="oc_ref") or None
        dg = solve_network(edges_with_ddg, reference=ref)
        if dg:
            dg_df = (pd.DataFrame(
                [{"ligand": k, "DeltaG (kcal/mol)": round(v[0], 2),
                  "uncertainty": round(v[1], 2)} for k, v in dg.items()])
                .sort_values("DeltaG (kcal/mol)").reset_index(drop=True))
            dg_df.insert(0, "rank", range(1, len(dg_df) + 1))
            st.dataframe(dg_df, hide_index=True, width="stretch")
            st.download_button("Download per-ligand DeltaG CSV",
                               dg_df.to_csv(index=False),
                               file_name=f"{wd.name}_dG.csv", mime="text/csv")
            st.markdown("#### 4) Cycle-closure QC")
            res = cycle_closure_residuals(edges_with_ddg, dg)
            if res:
                qc = pd.DataFrame(
                    [{"cycle": " -> ".join(c), "residual (kcal/mol)": round(r, 2),
                      "flag": "OK" if abs(r) <= 1.0 else "CHECK (>1)"}
                     for c, r in res])
                st.dataframe(qc, hide_index=True, width="stretch")
                worst = max(abs(r) for _, r in res)
                if worst > 1.0:
                    st.warning(f"Largest cycle residual {worst:.2f} kcal/mol "
                               "exceeds 1.0 - inspect those edges.")
                else:
                    st.success("All cycle residuals within 1.0 kcal/mol.")
            else:
                st.caption("No closed cycles in this network.")
    else:
        st.caption("Need an edge/DDG report to solve the network.")

with tab_script:
    st.code(oc.preview_script(layout, settings), language="bash")
