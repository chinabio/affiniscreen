"""
6_Results_Compare.py  --  Results: Compare / Rank (v2.5.0, Phase 3).

Multi-molecule view: rank a batch/campaign of ligands by predicted binding,
compare ΔG/ΔΔG, chart them, and (optionally) score against experimental
values pasted/uploaded by the user.

Sources per ligand subdir (auto-detected via results_lib.list_ligand_subdirs):
  * MM-GBSA  -> mmgbsa_status()  -> ΔG_bind
  * FEP/ABFE -> fep_headline()   -> ΔG_bind / ΔΔG
  * OpenFE   -> parse_openfe_result() (*_result.json) -> ΔG_bind
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="Results: Compare", layout="wide",
                   page_icon="trophy")

from amber_md.gui.common import dir_picker
from amber_md.gui import results_lib as rl
from amber_md.gui import routing_help as rh
from amber_md.gui import amber_campaign as oc_amber
from amber_md.gui import openfe_campaign as oc_openfe
from amber_md.gui.fep_common import solve_network, cycle_closure_residuals


def _render_openfe_report_section(parent):
    """Render the OpenFE report UI. Works for any dir with *_result.json
    or <ligand>_rep*/results.json, independent of list_ligand_subdirs."""
    # ---------------------------------------------------------------------------
    # OpenFE FEP/ABFE HTML report generator (v2.5.0)
    #   Auto-discovers *_result.json per ligand (ANY replicate count), writes one
    #   self-contained HTML per ligand + a campaign index. Cross-check vs an
    #   external reference (e.g. FEP+) is shown ONLY when explicitly supplied.
    # ---------------------------------------------------------------------------
    st.divider()
    st.subheader("Generate OpenFE FEP/ABFE report")
    st.caption("Scan the campaign directory for *_result.json (any number of "
               "replicates per ligand) and build a single-page HTML report per "
               "ligand plus a campaign summary.")

    with st.expander("Report options", expanded=False):
        rep_out = st.text_input(
            "Output directory (blank = <campaign>/_report)", "",
            key="ofe_report_out")
        st.markdown("**External reference (optional)** -- normally leave blank; "
                    "there is usually no FEP+ data. One per line as "
                    "`ligand_substring, label, value`, e.g. `12944901, FEP+, -7.56`.")
        ref_text = st.text_area("References", "", key="ofe_report_refs",
                                height=80, placeholder="12944901, FEP+, -7.56")

    if st.button("Build report", type="primary", key="ofe_report_go"):
        from amber_md.gui import openfe_report as ofr
        # parse optional references
        reference = {}
        for line in (ref_text or "").splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 3:
                try:
                    reference[parts[0]] = (parts[1], float(parts[2]))
                except ValueError:
                    st.warning(f"Ignored bad reference line: {line!r}")
        out_dir = rep_out.strip() or None
        try:
            with st.spinner("Discovering results and rendering HTML..."):
                res = ofr.build_campaign_report(parent, out_dir=out_dir,
                                                reference=reference or None)
            rows = res["rows"]
            if not rows:
                st.warning("No *_result.json files found under the campaign "
                           "directory.")
            else:
                mode = res.get("mode", "abfe")
                noun = "transformation(s)" if mode == "rbfe" else "ligand(s)"
                st.success(f"Report written for {len(rows)} {noun}. "
                           f"Index: {res['index']}")
                if mode == "rbfe":
                    # RBFE rows: ligA/ligB/ddg_pred/sem/sd/dG_complex/dG_solvent
                    summ = pd.DataFrame([{
                        "transformation": f"{r['ligA']} \u2192 {r['ligB']}",
                        "n reps": r.get("n_reps"),
                        "ΔΔG pred": (None if r.get("ddg_pred") is None
                                     else round(r["ddg_pred"], 2)),
                        "SEM": (None if r.get("sem") is None
                                else round(r["sem"], 2)),
                        "ΔG complex": (None if r.get("dG_complex") is None
                                       else round(r["dG_complex"], 2)),
                        "ΔG solvent": (None if r.get("dG_solvent") is None
                                       else round(r["dG_solvent"], 2)),
                        "ΔΔG exp": (None if r.get("ddg_exp") is None
                                    else round(r["ddg_exp"], 2)),
                        "|error|": (None if r.get("abs_err") is None
                                    else round(r["abs_err"], 2)),
                    } for r in rows])
                    st.dataframe(summ, hide_index=True, width="stretch")
                    stats = res.get("stats") or {}
                    if stats:
                        _r = stats.get("pearson_r")
                        _t = stats.get("kendall_tau")
                        st.caption(
                            "vs experiment (n={n}): MUE={mue:.2f}, RMSE={rmse:.2f}, "
                            "Pearson R={pr}, Kendall tau={kt}".format(
                                n=stats.get("n"),
                                mue=stats.get("mue", float("nan")),
                                rmse=stats.get("rmse", float("nan")),
                                pr=("%.2f" % _r) if _r is not None else "n/a",
                                kt=("%.2f" % _t) if _t is not None else "n/a"))
                else:
                    summ = pd.DataFrame([{
                        "ligand": r["key"], "n reps": r["stats"]["n"],
                        "mean ΔG": (None if r["stats"]["mean"] is None
                                    else round(r["stats"]["mean"], 2)),
                        "SEM": (None if r["stats"]["sem"] is None
                                else round(r["stats"]["sem"], 2)),
                        "SD": (None if r["stats"]["sd"] is None
                               else round(r["stats"]["sd"], 2)),
                        "status": r["status"],
                    } for r in rows])
                    st.dataframe(summ, hide_index=True, width="stretch")
                # offer the campaign index for download
                try:
                    idx_html = Path(res["index"]).read_text()
                    st.download_button("Download campaign index.html", idx_html,
                                       file_name="index.html", mime="text/html",
                                       key="ofe_report_dl")
                except Exception:
                    pass
                st.info(f"Per-ligand HTML files are in: "
                        f"{Path(res['index']).parent}")
        except Exception as e:  # noqa: BLE001
            st.error(f"Report generation failed: {e}")



st.title("Results - Compare & Rank")
st.caption("Rank all ligands in a campaign and score them against experimental "
           "values. Works for MM-GBSA and OpenFE ABFE/RBFE results.")
st.caption("Rank a campaign of ligands by predicted binding and compare them. "
           "Lower ΔG = stronger predicted binder.")

exp = st.session_state.get("experiment", {})
default_wd = exp.get("work_dir", "~/Run_dir")
with st.sidebar:
    st.header("Campaign directory")
    parent = dir_picker("Batch / parent directory", "rc_dir",
                        default_dir=default_wd)
parent = Path(parent).expanduser() if parent else None

if not parent or not parent.exists():
    st.info("Select a parent directory containing per-ligand subdirectories.")
    st.stop()

lig_dirs = rl.list_ligand_subdirs(parent)
if not lig_dirs:
    st.warning(f"No ranked ligand subdirectories found under `{parent}` "
               "(expected lig_*/ or subdirs with fep/ or mmgbsa/). "
               "The OpenFE report generator below still works on *_result.json "
               "or <ligand>_rep*/results.json layouts.")
    _render_openfe_report_section(parent)
    st.stop()

st.subheader(f"Campaign: `{parent.name}`  ({len(lig_dirs)} ligands)")

# final62: routing help (which page does what) -- collapsed by default.
rh.render_routing_help(expanded=False)

# final62: RBFE NETWORK rollup. Per-ligand subdir scanning cannot rank a
# perturbation NETWORK (per-edge ddG must be SOLVED into per-ligand dG). Detect
# a network here and surface the solved ranking + cycle-closure QC inline, for
# BOTH engines, using the same solver as the FEP Campaign page.
def _network_edges_with_ddg(parent: Path):
    """Return ([(a,b,ddG,err)], engine_label) if parent is an RBFE network."""
    amb = oc_amber.CampaignLayout.at(parent)
    amb_edges = oc_amber.edge_ddg_table(amb)
    if amb_edges:
        return amb_edges, "Amber"
    ofe = oc_openfe.CampaignLayout.at(parent)
    tsv = ofe.dg_tsv
    if tsv.exists():
        try:
            import csv as _csv
            rows = []
            with open(tsv) as fh:
                rd = _csv.DictReader(fh, delimiter="\t")
                for row in rd:
                    keys = {k.lower(): k for k in row}
                    a = row.get(keys.get("ligand_a", ""), "")
                    b = row.get(keys.get("ligand_b", ""), "")
                    dgk = next((keys[k] for k in keys if "ddg" in k), None)
                    erk = next((keys[k] for k in keys
                                if "unc" in k or "err" in k), None)
                    if a and b and dgk:
                        rows.append((a, b, float(row[dgk]),
                                     float(row[erk]) if erk and row.get(erk)
                                     else 0.5))
            if rows:
                return rows, "OpenFE"
        except Exception:
            pass
    return [], ""


_net_edges, _net_engine = _network_edges_with_ddg(parent)
if _net_edges:
    st.markdown("### RBFE network (solved -> per-ligand ΔG)")
    st.caption(f"Detected a {_net_engine} perturbation network "
               f"({len(_net_edges)} edge(s) with ΔΔG). Solving into per-ligand "
               "ΔG with cycle-closure QC (same math as the FEP Campaign page).")
    _ref = st.text_input("Reference ligand (pinned to 0; blank = first)", "",
                         key="rc_net_ref") or None
    _dg = solve_network(_net_edges, reference=_ref)
    if _dg:
        _dg_df = (pd.DataFrame(
            [{"ligand": k, "ΔG (kcal/mol)": round(v[0], 2),
              "± err": round(v[1], 2), "method": "RBFE",
              "engine": _net_engine} for k, v in _dg.items()])
            .sort_values("ΔG (kcal/mol)").reset_index(drop=True))
        _dg_df.insert(0, "rank", range(1, len(_dg_df) + 1))
        st.dataframe(_dg_df, hide_index=True, width="stretch")
        st.download_button("Download network per-ligand ΔG CSV",
                           _dg_df.to_csv(index=False),
                           file_name=f"{parent.name}_rbfe_network_dG.csv",
                           mime="text/csv")
        _res = cycle_closure_residuals(_net_edges, _dg)
        if _res:
            _worst = max(abs(r) for _, r in _res)
            _qc = pd.DataFrame(
                [{"cycle": " -> ".join(c), "residual (kcal/mol)": round(r, 2),
                  "flag": "OK" if abs(r) <= 1.0 else "CHECK (>1)"}
                 for c, r in _res])
            with st.expander(f"Cycle-closure QC (worst {_worst:.2f} kcal/mol)"):
                st.dataframe(_qc, hide_index=True, width="stretch")
    st.markdown("---")


# ---- collect one row per ligand (best available result) ----
# final62: every row now carries explicit ENGINE + refined METHOD so a mixed
# campaign (Amber + OpenFE, MM-GBSA + ABFE) is unambiguous in the table.
def _method_engine(kind: str) -> tuple[str, str]:
    """Map a fep_headline 'kind' string to (method, engine)."""
    k = kind.lower()
    if "rbfe" in k or "relative" in k:
        return "RBFE", "Amber"
    if "abfe" in k or "absolute" in k or "dg_bind" in k:
        return "ABFE", "Amber"
    if "complex leg" in k or "solvent leg" in k:
        return "FEP (partial)", "Amber"
    return "FEP", "Amber"


def _collect(wd: Path) -> dict:
    name = wd.name
    # OpenFE first (richest) -- *_result.json => OpenFE engine.
    ofe = list(wd.rglob("*_result.json")) + list(wd.rglob("results.json"))
    if ofe:
        r = rl.parse_openfe_result(ofe[0])
        if r and r.get("estimate_kcal_mol") is not None:
            unc = r.get("uncertainty_kcal_mol")
            if not unc:
                errs = [v.get("unit_estimate_error") or 0.0
                        for v in r["legs"].values()]
                unc = math.sqrt(sum(e * e for e in errs)) if errs else None
            ofe_method = "RBFE" if "rbfe" in ofe[0].name.lower() else "ABFE"
            return {"ligand": name, "ΔG (kcal/mol)": r["estimate_kcal_mol"],
                    "± err": unc, "method": ofe_method, "engine": "OpenFE",
                    "status": "DONE"}
    # Amber FEP (fep/ABFE_RESULT.json or per-leg summaries)
    fep = wd / "fep"
    if fep.exists():
        val, kind, _ = rl.fep_headline(fep)
        if val is not None:
            method, engine = _method_engine(kind)
            return {"ligand": name, "ΔG (kcal/mol)": val, "± err": None,
                    "method": method, "engine": engine, "status": "DONE",
                    "detail": kind}
    # MM-GBSA (engine via the locked detector: mmgbsa/engine.json marker,
    # else OpenMM-exclusive MD artifacts, else Amber) -- see results_lib.
    status, dg, full = rl.mmgbsa_status(wd)
    engine = rl.mmgbsa_engine(wd)
    return {"ligand": name, "ΔG (kcal/mol)": dg,
            "± err": (full.get("sem") if full else None),
            "method": "MM-GBSA", "engine": engine, "status": status}


rows = [_collect(d) for d in lig_dirs]
df = pd.DataFrame(rows)
done = df[df["ΔG (kcal/mol)"].notna()].copy()

# ---- summary ----
m1, m2, m3 = st.columns(3)
m1.metric("Ligands", len(df))
m2.metric("With ΔG", len(done))
m3.metric("Pending", len(df) - len(done))

# ============================ EXPERIMENTAL (optional) ============================
with st.expander("Add experimental ΔG (optional, for RMSE/ranking score)"):
    st.caption("Paste CSV with columns: ligand,exp_dG  (kcal/mol).")
    exp_text = st.text_area("ligand,exp_dG", "", key="rc_exp") or ""
    exp_map = {}
    if exp_text.strip():
        for line in exp_text.strip().splitlines():
            parts = [p.strip() for p in line.replace("\t", ",").split(",")]
            if len(parts) >= 2 and parts[0].lower() != "ligand":
                try:
                    exp_map[parts[0]] = float(parts[1])
                except ValueError:
                    pass
        st.caption(f"Parsed {len(exp_map)} experimental value(s).")

if not len(done):
    st.info("No finished ΔG values yet. The full status table is below.")
    st.dataframe(df, hide_index=True, width="stretch")
    st.stop()

# attach experimental + signed error
if exp_map:
    done["exp ΔG"] = done["ligand"].map(exp_map)
    done["|Δ|"] = (done["ΔG (kcal/mol)"] - done["exp ΔG"]).abs()

# ---- filters / sort ----
f1, f2 = st.columns(2)
sort_asc = f1.checkbox("Sort ascending (strongest first)", value=True,
                       key="rc_sortasc")
max_err = f2.number_input("Hide |Δ| above (kcal/mol, 0 = no filter)",
                          0.0, 50.0, 0.0, 0.5, key="rc_maxerr")

view = done.sort_values("ΔG (kcal/mol)", ascending=sort_asc).reset_index(drop=True)
if exp_map and max_err > 0 and "|Δ|" in view.columns:
    view = view[view["|Δ|"].fillna(0) <= max_err]
view.insert(0, "rank", range(1, len(view) + 1))

st.markdown("### Ranking")
st.dataframe(view, hide_index=True, width="stretch")

# ---- RMSE vs experiment ----
if exp_map and "|Δ|" in done.columns:
    paired = done.dropna(subset=["exp ΔG"])
    if len(paired) >= 2:
        rmse = math.sqrt((paired["|Δ|"] ** 2).mean())
        mae = paired["|Δ|"].mean()
        try:
            from scipy.stats import pearsonr
            r_val = pearsonr(paired["ΔG (kcal/mol)"], paired["exp ΔG"])[0]
        except Exception:
            r_val = paired["ΔG (kcal/mol)"].corr(paired["exp ΔG"])
        c1, c2, c3 = st.columns(3)
        c1.metric("RMSE vs exp", f"{rmse:.2f}")
        c2.metric("MAE", f"{mae:.2f}")
        c3.metric("Pearson r", f"{r_val:.2f}")

# ---- bar chart ----
if len(view) >= 2:
    try:
        import altair as alt
        chart = alt.Chart(view).mark_bar().encode(
            x=alt.X("ΔG (kcal/mol):Q"),
            y=alt.Y("ligand:N", sort="x", title=None),
            color=alt.Color("ΔG (kcal/mol):Q",
                            scale=alt.Scale(scheme="redblue", reverse=True)),
            tooltip=list(view.columns),
        ).properties(height=max(200, 28 * len(view)))
        st.altair_chart(chart, width="stretch")
    except Exception as e:               # noqa: BLE001
        st.caption(f"(chart unavailable: {e})")

# ---- scatter vs experiment ----
if exp_map and "exp ΔG" in view.columns and view["exp ΔG"].notna().sum() >= 2:
    try:
        import altair as alt
        sc = view.dropna(subset=["exp ΔG"])
        base = alt.Chart(sc).encode(
            x=alt.X("exp ΔG:Q", title="Experimental ΔG"),
            y=alt.Y("ΔG (kcal/mol):Q", title="Predicted ΔG"),
            tooltip=["ligand", "ΔG (kcal/mol)", "exp ΔG", "|Δ|"])
        pts = base.mark_circle(size=90)
        lo = float(min(sc["exp ΔG"].min(), sc["ΔG (kcal/mol)"].min())) - 1
        hi = float(max(sc["exp ΔG"].max(), sc["ΔG (kcal/mol)"].max())) + 1
        diag = alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo, hi]})).mark_line(
            strokeDash=[4, 4], color="gray").encode(x="x", y="y")
        st.altair_chart((diag + pts).properties(height=360), width="stretch")
    except Exception:
        pass

# ---- export ----
st.download_button("Download ranking CSV", view.to_csv(index=False),
                   file_name=f"{parent.name}_ranking.csv", mime="text/csv")

st.markdown("### Full status (all ligands)")
st.dataframe(df, hide_index=True, width="stretch")


# NOTE: the 'Promote MM-GBSA -> FEP' action was removed from the GUI (v2.6.0)
# because it scaffolded and submitted **Amber** ABFE/RBFE runs via
# amber_md.fep_driver -- a workflow no longer exposed in the GUI. The helper
# module amber_md.gui.promote_fep is retained for programmatic / CLI use.

# Render the OpenFE report section on the normal path too.
_render_openfe_report_section(parent)
