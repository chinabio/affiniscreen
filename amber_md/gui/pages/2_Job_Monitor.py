"""
2_Job_Monitor.py  --  Live job monitor (v2.5.0, Phase 2).

Aggregates jobs from THREE sources and presents one unified, auto-refreshing
table with status, progress, and a log preview:

  1. LSF queue        -> get_lsf_jobs()  (PEND/RUN/DONE/EXIT)
  2. Detached launches -> st.session_state["launches"] (login-node PIDs + logs
     recorded by the Setup & Launch wizard)
  3. OpenFE edges     -> per-edge result JSON status via edge_repeat_status()
     for any work_dir that contains a network_setup/transformations layout.

Nothing here changes engine behaviour; it is read-only monitoring built on the
existing helpers (get_lsf_jobs, tail_file, is_pid_alive, parse_md_progress,
edge_repeat_status, list_transformations).
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import time
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="Job Monitor", layout="wide", page_icon="bar_chart")

from amber_md.gui.common import (
    get_lsf_jobs, tail_file, is_pid_alive, parse_md_progress,
)
from amber_md.gui.openfe_common import (
    list_transformations, edge_repeat_status,
)

# ----------------------------------------------------------------------------
# Status model: normalise every source to a small common vocabulary.
# ----------------------------------------------------------------------------
STATUS_ICON = {
    "queued":    "QUEUED",
    "running":   "RUNNING",
    "completed": "DONE",
    "failed":    "FAILED",
    "unknown":   "?",
}
# colour hints for the dataframe styling
STATUS_COLOR = {
    "queued":   "#fff3cd",
    "running":  "#d7e9f7",
    "completed": "#cfe8cf",
    "failed":   "#f7d7d7",
    "unknown":  "#eeeeee",
}

_LSF_MAP = {
    "PEND": "queued", "WAIT": "queued", "PROV": "queued",
    "RUN": "running", "DONE": "completed", "EXIT": "failed",
    "PSUSP": "queued", "USUSP": "queued", "SSUSP": "queued",
    "UNKWN": "unknown", "ZOMBI": "failed",
}


def _norm_lsf(stat: str) -> str:
    return _LSF_MAP.get((stat or "").upper(), "unknown")


# ----------------------------------------------------------------------------
# Launch registry: persist launches recorded by the wizard so the monitor can
# track login-node detached processes (those never appear in bjobs).
# The wizard writes st.session_state["last_launch"]; we accumulate them into a
# durable list so multiple launches are all tracked.
# ----------------------------------------------------------------------------
def _ingest_last_launch():
    last = st.session_state.get("last_launch")
    if not last:
        return
    reg = st.session_state.setdefault("launches", [])
    # dedupe by (work_dir, ts)
    sig = (last.get("work_dir"), last.get("ts"))
    if not any((r.get("work_dir"), r.get("ts")) == sig for r in reg):
        reg.append(dict(last))


def _detached_rows() -> list[dict]:
    rows = []
    for r in st.session_state.get("launches", []):
        method = r.get("method", "?")
        engine = r.get("engine", "?")
        wd = r.get("work_dir", "")
        pids = r.get("pids", []) or []
        logs = r.get("logs", []) or []
        for i, pid in enumerate(pids):
            log = logs[i] if i < len(logs) else ""
            alive = is_pid_alive(pid)
            tail = tail_file(log, 60) if log else ""
            prog = ""
            # Status precedence:
            #   1. PID still alive            -> running (authoritative)
            #   2. error/traceback in log     -> failed
            #   3. explicit success marker    -> completed
            #   4. process gone, no markers   -> completed (best guess)
            if alive:
                status = "running"
            elif "Traceback" in tail or "Error:" in tail or "ERROR" in tail:
                status = "failed"
            elif ("Duration:" in tail or "Done" in tail or "FINAL" in tail
                  or "complete" in tail.lower()):
                status = "completed"
            else:
                status = "completed"
            rows.append({
                "source": "login-PID",
                "method": method, "engine": engine,
                "name": f"{method} launch pid={pid}",
                "status": status, "progress": prog,
                "work_dir": wd, "log": log, "id": str(pid),
            })
    return rows


def _lsf_rows() -> list[dict]:
    rows = []
    for j in get_lsf_jobs():
        rows.append({
            "source": "LSF",
            "method": "", "engine": "",
            "name": j.get("name", ""),
            "status": _norm_lsf(j.get("stat", "")),
            "progress": "",
            "work_dir": "", "log": "", "id": j.get("jobid", ""),
            "queue": j.get("queue", ""), "submitted": j.get("submitted", ""),
        })
    return rows


def _openfe_edge_rows() -> list[dict]:
    """For each registered OpenFE work_dir, report per-edge result status."""
    rows = []
    seen = set()
    for r in st.session_state.get("launches", []):
        if "OpenFE" not in r.get("engine", ""):
            continue
        wd = Path(r.get("work_dir", "")).expanduser()
        # final60: ABFE plans into abfe_setup/, RBFE into network_setup/.
        # Use the same resolver the campaign submitter uses so ABFE edges show
        # up in the Monitor instead of silently reporting nothing.
        from amber_md.gui.openfe_campaign import resolve_network_dir
        net = resolve_network_dir(wd)
        if not net.is_dir() or str(net) in seen:
            continue
        seen.add(str(net))
        edges = list_transformations(net)
        for e in edges:
            rjson = e.with_name(e.stem + "_result.json")
            status = edge_repeat_status(rjson, e.parent)
            norm = {"DONE": "completed", "running": "running",
                    "analyzing": "running", "not started": "queued"}.get(
                        status, "unknown")
            rows.append({
                "source": "OpenFE-edge",
                "method": r.get("method", "RBFE"),
                "engine": r.get("engine", ""),
                "name": e.stem,
                "status": norm, "progress": "",
                "work_dir": str(wd), "log": "", "id": e.stem,
            })
    return rows


def _amber_progress(work_dir: str) -> str:
    """Best-effort Amber MD progress from a prod.out under the work dir."""
    wd = Path(work_dir).expanduser()
    if not wd.exists():
        return ""
    for cand in wd.rglob("*prod*.out"):
        pr = parse_md_progress(cand)
        if pr:
            step, ps = pr
            return f"{ps:.0f} ps"
    return ""


# ============================ Header & controls ============================
st.title("Job Monitor")
st.caption("See what's running, queued, or finished across all your jobs -- "
           "and tail their logs.")

# --- Disk-scan completion view (session-independent) -------------------------
# Reads finished jobs straight from a work_dir on disk, so completed runs are
# visible even after the app restarts / a new browser session, and for jobs
# launched from the CLI. Validates OpenFE result JSONs (tiny/failed files are
# reported as FAILED, not completed).
with st.expander("\U0001F50D Scan a directory for finished jobs "
                 "(works across restarts / CLI launches)", expanded=False):
    from amber_md.gui.openfe_common import scan_results_dir
    from amber_md.gui.common import dir_picker
    _exp = st.session_state.get("experiment", {})
    _scan_dir = dir_picker("Work / results directory to scan", "mon_scan_dir",
                           default_dir=_exp.get("work_dir", "~/Run_dir"))
    if st.button("Scan now", key="mon_scan_btn") and _scan_dir:
        res = scan_results_dir(_scan_dir)
        if res.get("error"):
            st.error(res["error"])
        else:
            c = res["counts"]
            cols = st.columns(4)
            cols[0].metric("Completed", c.get("completed", 0))
            cols[1].metric("Failed", c.get("failed", 0))
            cols[2].metric("Running", c.get("running", 0))
            cols[3].metric("Unreadable", c.get("unreadable", 0))
            import pandas as _pd
            _icon = {"completed": "\u2705", "failed": "\u274c",
                     "running": "\U0001F535", "unreadable": "\u2754"}
            df = _pd.DataFrame([{
                "": _icon.get(r["status"], "?"),
                "status": r["status"], "kind": r["kind"],
                "name": r["name"], "detail": r["detail"], "path": r["path"],
            } for r in res["rows"]])
            if df.empty:
                st.info("No result.json or MM-GBSA .dat files found under "
                        f"`{_scan_dir}`.")
            else:
                st.dataframe(df, hide_index=True, width="stretch")
                done = df[df["status"] == "completed"]
                st.download_button(
                    "Download completed list (CSV)",
                    done.to_csv(index=False),
                    file_name="completed_jobs.csv", mime="text/csv",
                    disabled=done.empty, key="mon_scan_csv")
                if (df["status"] == "failed").any():
                    st.warning("Some result files exist but look **failed** "
                               "(no estimate / carries an error). They are NOT "
                               "counted as completed -- inspect the path.")
                    # Show the REAL failure reason mined from each failed result.
                    for _r in res["rows"]:
                        if _r["status"] != "failed":
                            continue
                        _f = _r.get("failure")
                        with st.expander(
                                f"\u274c Why did **{_r['name']}** fail?",
                                expanded=False):
                            if not _f:
                                st.write(_r["detail"]); 
                                st.caption(_r["path"]); 
                                continue
                            st.markdown(
                                f"**{_f.get('headline','(no headline)')}**")
                            st.caption(
                                f"{_f.get('n_failed',0)}/{_f.get('n_units',0)} "
                                f"protocol unit(s) failed, "
                                f"{_f.get('n_ok',0)} ok.")
                            for _u in _f.get("failures", []):
                                st.markdown(
                                    f"- `{_u.get('exc_type','')}` "
                                    f"{_u.get('exc_msg','')}  \n"
                                    f"  *unit:* {_u.get('unit','')}  \n"
                                    f"  *origin:* `{_u.get('origin','') or 'n/a'}`")
                            _tb = ""
                            if _f.get("failures"):
                                _tb = _f["failures"][0].get("traceback", "")
                            _tb = _tb or _f.get("top_traceback", "")
                            if _tb:
                                with st.expander("Show full traceback",
                                                 expanded=False):
                                    st.code(_tb, language="text")
                            st.caption(_r["path"])


st.caption("Live status across scheduler (LSF) jobs, detached login-node launches, and OpenFE "
           "edges. Expand a row for the log tail.")

_ingest_last_launch()

ctop1, ctop2, ctop3, ctop4 = st.columns([1.2, 1, 1, 1])
auto = ctop1.checkbox("Auto-refresh", value=True, key="mon_auto")
interval = ctop2.number_input("Every (s)", 3, 60, 10, key="mon_int")
ctop3.button("Refresh now", key="mon_refresh")   # any click triggers a rerun
if ctop4.button("Clear finished from registry", key="mon_clear"):
    st.session_state["launches"] = [
        r for r in st.session_state.get("launches", [])
        if any(is_pid_alive(p) for p in (r.get("pids") or []))
    ]

# filters
f1, f2, f3 = st.columns(3)
src_filter = f1.multiselect(
    "Source", ["LSF", "login-PID", "OpenFE-edge"],
    default=["LSF", "login-PID", "OpenFE-edge"], key="mon_src")
status_filter = f2.multiselect(
    "Status", list(STATUS_ICON.keys()),
    default=["queued", "running", "failed", "completed"], key="mon_stat")
search = f3.text_input("Search name", "", key="mon_search")

# ============================ Gather rows ============================
rows = []
if "LSF" in src_filter:
    rows += _lsf_rows()
if "login-PID" in src_filter:
    rows += _detached_rows()
if "OpenFE-edge" in src_filter:
    rows += _openfe_edge_rows()

# enrich Amber login-PID rows with MD progress
for r in rows:
    if r["source"] == "login-PID" and not r["progress"] and r["work_dir"]:
        r["progress"] = _amber_progress(r["work_dir"])

# apply filters
def _keep(r):
    if r["status"] not in status_filter:
        return False
    if search and search.lower() not in r["name"].lower():
        return False
    return True

rows = [r for r in rows if _keep(r)]

# ============================ Summary metrics ============================
counts = {k: 0 for k in STATUS_ICON}
for r in rows:
    counts[r["status"]] = counts.get(r["status"], 0) + 1
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Total", len(rows))
m2.metric("Running", counts["running"])
m3.metric("Queued", counts["queued"])
m4.metric("Completed", counts["completed"])
m5.metric("Failed", counts["failed"])

# ============================ Table ============================
if rows:
    df = pd.DataFrame([{
        "Status": STATUS_ICON.get(r["status"], "?"),
        "Source": r["source"],
        "Name": r["name"],
        "Method": r.get("method", ""),
        "Engine": r.get("engine", ""),
        "Progress": r.get("progress", ""),
        "ID": r.get("id", ""),
        "Queue": r.get("queue", ""),
    } for r in rows])

    def _row_style(row):
        # map the icon back to a status colour
        inv = {v: k for k, v in STATUS_ICON.items()}
        c = STATUS_COLOR.get(inv.get(row["Status"], "unknown"), "#fff")
        return [f"background-color: {c}"] * len(row)

    try:
        st.dataframe(df.style.apply(_row_style, axis=1),
                     hide_index=True, width="stretch")
    except Exception:
        st.dataframe(df, hide_index=True, width="stretch")

    # ---- per-job log preview ----
    st.subheader("Log preview")
    logged = [r for r in rows if r.get("log")]
    if logged:
        labels = [f'{r["name"]}  ({r["status"]})' for r in logged]
        pick = st.selectbox("Pick a job with a log", range(len(logged)),
                            format_func=lambda i: labels[i], key="mon_pick")
        chosen = logged[pick]
        st.caption(f'Log: `{chosen["log"]}`')
        st.code(tail_file(chosen["log"], 60) or "(empty)")
    else:
        st.caption("No login-node logs available for the current rows. "
                   "(LSF jobs write to their own .out/.err files.)")
else:
    st.info("No jobs match the current filters. Launch something from "
            "**Setup & Launch**, or widen the filters above.")

# ============================ Auto-refresh loop ============================
# Only rerun while there is at least one active (queued/running) job, to avoid
# spinning forever once everything is finished.
active = any(r["status"] in ("queued", "running") for r in rows)
if auto and active:
    time.sleep(int(interval))
    st.rerun()
elif auto and not active and rows:
    st.caption("All tracked jobs are finished -- auto-refresh paused.")
