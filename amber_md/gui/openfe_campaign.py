"""
openfe_campaign.py  --  Streamlit-free OpenFE campaign orchestration
(v2.5.0, Phase 4).

Extracts the proven Plan -> Run -> Gather -> Solve lifecycle from the legacy
7_RBFE_OpenFE.py Run/Results tabs into reusable, testable functions so the
unified Setup/Monitor flow can drive an OpenFE campaign end-to-end without the
legacy page.

Pure logic only (no `import streamlit`). Side effects (writing scripts,
submitting bsub, spawning gather) are isolated in small functions a UI layer
calls; everything else is queryable status.
"""
# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations

import time
import shlex
from dataclasses import dataclass
from pathlib import Path

from amber_md.gui.openfe_common import (
    OpenFESettings, make_edge_bsub_script, submit_bsub, gather_cmd,
    count_my_jobs, list_transformations, edge_repeat_status, MAX_GPU,
)

JOB_PREFIX = "ofe_camp"


# ---------------------------------------------------------------------------
# Standard directory layout for a campaign rooted at work_dir.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# final60 FIX: ABFE/OpenFE plans transformations into <wd>/abfe_setup/, while
# RBFE/OpenFE plans into <wd>/network_setup/. The campaign submitter, status
# snapshot and mapping-artifact finder all hard-coded "network_setup", so an
# ABFE campaign saw ZERO transformations and submitted nothing for ANY
# molecule. resolve_network_dir() picks whichever planning dir actually holds
# transformations, falling back to network_setup for backward compatibility.
# ---------------------------------------------------------------------------
_NETWORK_DIR_CANDIDATES = ("network_setup", "abfe_setup")


def resolve_network_dir(work_dir: Path) -> Path:
    """Return the planning dir under work_dir that contains transformations.

    Prefers a candidate whose transformations/ holds >=1 *.json (or, failing
    that, the bare dir holds *.json). If none qualifies, returns
    <work_dir>/network_setup so existing RBFE behaviour is preserved and the
    page can still show its "no transformations" hint.
    """
    wd = Path(work_dir).expanduser()
    for name in _NETWORK_DIR_CANDIDATES:
        cand = wd / name
        tdir = cand / "transformations"
        try:
            if tdir.is_dir() and any(tdir.glob("*.json")):
                return cand
            if cand.is_dir() and any(cand.glob("*.json")):
                return cand
        except OSError:
            continue
    return wd / "network_setup"


@dataclass
class CampaignLayout:
    work_dir: Path
    network_dir: Path
    results_dir: Path
    log_dir: Path
    dg_tsv: Path

    @classmethod
    def at(cls, work_dir: Path) -> "CampaignLayout":
        wd = Path(work_dir).expanduser()
        return cls(
            work_dir=wd,
            network_dir=resolve_network_dir(wd),
            results_dir=wd / "results",
            log_dir=wd / "logs",
            dg_tsv=wd / "rbfe_dg.tsv",
        )


def edges(layout: CampaignLayout) -> list[Path]:
    return list_transformations(layout.network_dir)


def task_list(layout: CampaignLayout, n_replicates: int):
    """All (edge_path, repeat_index) tasks for the campaign."""
    return [(e, r) for e in edges(layout) for r in range(n_replicates)]


def task_paths(layout: CampaignLayout, edge: Path, rep: int):
    name = f"{edge.stem}_rep{rep}"
    return name, (layout.results_dir / f"{name}_result.json",
                  layout.results_dir / f"{name}_work")


# ---------------------------------------------------------------------------
# Status snapshot (read-only) -- powers the Monitor and progress bars.
# ---------------------------------------------------------------------------
def campaign_status(layout: CampaignLayout, n_replicates: int) -> dict:
    rows, buckets = [], {"DONE": 0, "analyzing": 0, "running": 0,
                         "not started": 0}
    for edge, rep in task_list(layout, n_replicates):
        name, (res_json, ework) = task_paths(layout, edge, rep)
        stt = edge_repeat_status(res_json, ework)
        buckets[stt] = buckets.get(stt, 0) + 1
        rows.append({"job": name, "status": stt,
                     "result_json": str(res_json) if res_json.exists() else "",
                     "edge": edge.stem, "repeat": rep})
    total = len(rows)
    return {"rows": rows, "buckets": buckets, "total": total,
            "done": buckets.get("DONE", 0),
            "pct": (buckets.get("DONE", 0) / total) if total else 0.0}


# ---------------------------------------------------------------------------
# Submit edges, throttled. `sleep_fn` and `slots_fn` are injectable for tests.
# Returns (submitted, skipped, errors[list[str]]).
# ---------------------------------------------------------------------------
def submit_campaign(layout: CampaignLayout, s: OpenFESettings,
                    n_replicates: int,
                    progress_cb=None,
                    slots_fn=None,         # accepted for API compat; unused
                    submit_fn=None,
                    sleep_fn=time.sleep,   # accepted for API compat; unused
                    poll_interval: int = 20,
                    throttle_n: int = 0,   # 0/None = OFF (default, unchanged)
                    throttle_poll: int = 30):
    """Submit every (edge, repeat) job.

    Default (throttle_n falsy): submit everything at once -- the LSF queue
    holds pending jobs. This is the v2.5.0 behaviour and is UNCHANGED.

    Optional (throttle_n >= 1): instead of submitting directly, write all the
    per-edge #BSUB scripts plus a joblist, then launch the detached login-node
    throttler (amber_md.throttle_submit) under nohup. The throttler keeps no
    more than throttle_n of THIS campaign's jobs in the queue at a time and
    SURVIVES this page being closed. Returns (n_planned, skipped, errors) with
    a marker error-list entry describing the handoff.
    """
    layout.results_dir.mkdir(parents=True, exist_ok=True)
    layout.log_dir.mkdir(parents=True, exist_ok=True)
    slots_fn = slots_fn or (lambda: count_my_jobs(JOB_PREFIX))
    submit_fn = submit_fn or submit_bsub

    tasks = task_list(layout, n_replicates)

    # ----- optional throttled handoff -------------------------------------
    if throttle_n and int(throttle_n) >= 1:
        return _submit_campaign_throttled(
            layout, s, tasks, int(throttle_n), int(throttle_poll), progress_cb)

    # ----- default: submit everything at once (UNCHANGED) -----------------
    submitted, skipped, errors = 0, 0, []
    for i, (edge, rep) in enumerate(tasks):
        name, (res_json, ework) = task_paths(layout, edge, rep)
        if res_json.exists():
            skipped += 1
            if progress_cb:
                progress_cb(i + 1, len(tasks), submitted, skipped)
            continue
        job_name = f"{JOB_PREFIX}_{name}"[:60]
        script = make_edge_bsub_script(edge, res_json, ework, s, job_name,
                                       layout.log_dir)
        script_path = layout.log_dir / f"submit_{name}.sh"
        script_path.write_text(script)
        ok, msg = submit_fn(script_path)
        if ok:
            submitted += 1
        else:
            errors.append(f"{name}: {msg}")
        if progress_cb:
            progress_cb(i + 1, len(tasks), submitted, skipped)
    return submitted, skipped, errors


def _submit_campaign_throttled(layout, s, tasks, throttle_n, throttle_poll,
                               progress_cb):
    """Write per-edge bsub scripts + a joblist, then launch the detached
    login-node throttler. Each joblist entry is `bsub < <script>` so the
    throttler's bjobs-based counter governs how many are in the queue.
    """
    import json as _json
    import sys as _sys
    import subprocess as _subprocess
    from datetime import datetime as _dt

    joblist = []
    skipped = 0
    for i, (edge, rep) in enumerate(tasks):
        name, (res_json, ework) = task_paths(layout, edge, rep)
        if res_json.exists():
            skipped += 1
            if progress_cb:
                progress_cb(i + 1, len(tasks), 0, skipped)
            continue
        job_name = f"{JOB_PREFIX}_{name}"[:60]
        script = make_edge_bsub_script(edge, res_json, ework, s, job_name,
                                       layout.log_dir)
        script_path = layout.log_dir / f"submit_{name}.sh"
        script_path.write_text(script)
        # `bsub < script` via a shell so the throttler submits exactly one job.
        joblist.append({
            "cmd": ["bash", "-lc", f"bsub < {shlex.quote(str(script_path))}"],
            "cwd": str(layout.log_dir),
            "log": str(layout.log_dir / f"throttled_{name}.log"),
        })
        if progress_cb:
            progress_cb(i + 1, len(tasks), len(joblist), skipped)

    if not joblist:
        return 0, skipped, ["All tasks already have results; nothing to submit."]

    joblist_path = layout.log_dir / "throttle_joblist.json"
    joblist_path.write_text(_json.dumps(joblist, indent=2))
    throttle_log = layout.log_dir / "throttle.log"

    cmd = [_sys.executable, "-m", "amber_md.throttle_submit",
           "--joblist", str(joblist_path),
           "--max-inflight", str(throttle_n),
           "--name-prefix", JOB_PREFIX,
           "--poll", str(throttle_poll),
           "--logfile", str(throttle_log)]
    # nohup + detached so it outlives the Streamlit page.
    try:
        with open(layout.log_dir / "throttle_nohup.out", "a") as _out:
            _out.write(f"\n# nohup throttler at {_dt.now().isoformat()}\n")
            _out.write(f"# {' '.join(cmd)}\n\n"); _out.flush()
            _subprocess.Popen(cmd, cwd=str(layout.log_dir),
                              stdout=_out, stderr=_subprocess.STDOUT,
                              start_new_session=True)
    except Exception as e:  # noqa: BLE001
        return 0, skipped, [f"Failed to launch throttler: {e}"]

    msg = (f"THROTTLED: handed {len(joblist)} job(s) to the login-node "
           f"throttler (max {throttle_n} in queue at a time). It runs under "
           f"nohup and survives closing this page. Progress: {throttle_log}")
    return len(joblist), skipped, [msg]


def preview_script(layout: CampaignLayout, s: OpenFESettings) -> str:
    es = edges(layout)
    if not es:
        return "# No transformations planned yet."
    sample = es[0]
    name = f"{sample.stem}_rep0"
    return make_edge_bsub_script(
        sample, layout.results_dir / f"{name}_result.json",
        layout.results_dir / f"{name}_work", s,
        f"{JOB_PREFIX}_{name}"[:60], layout.log_dir)


def find_mapping_artifacts(layout: "CampaignLayout") -> dict:
    """Locate atom-map inspection files under the work dir.

    Looks for rbfe_map.py outputs (edges.csv / nodes.csv / mapping.json /
    diagnostics.txt) in common locations, plus the OpenFE planned-network dir.
    Returns {kind: Path} for whatever exists.
    """
    wd = layout.work_dir
    found: dict = {}
    search_dirs = [wd, wd / "network_setup", wd / "abfe_setup", wd / "map", wd / "mapping",
                   layout.network_dir]
    for d in search_dirs:
        if not d.is_dir():
            continue
        for kind, fname in (("edges", "edges.csv"), ("nodes", "nodes.csv"),
                            ("mapping", "mapping.json"),
                            ("diagnostics", "diagnostics.txt"),
                            ("graphml", "map.graphml")):
            p = d / fname
            if kind not in found and p.exists():
                found[kind] = p
    if layout.network_dir.is_dir():
        found["transformations_dir"] = layout.network_dir
    return found


def read_edges_table(edges_csv: "Path") -> list:
    """Parse edges.csv into a list of dict rows (no pandas dependency)."""
    import csv
    rows = []
    with open(edges_csv, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def mapping_network_mermaid(edges_rows: list) -> str:
    """Build a Mermaid graph definition from edge rows (lig_a -- lig_b)."""
    lines = ["graph LR"]
    for e in edges_rows:
        a = (e.get("lig_a") or "").replace('"', "'")
        b = (e.get("lig_b") or "").replace('"', "'")
        if not a or not b:
            continue
        score = e.get("score") or e.get("cost") or ""
        label = f"|{score}|" if score else ""
        lines.append(f'    "{a}" ---{label} "{b}"')
    return "\n".join(lines)


def gather_command(layout: CampaignLayout, report: str = "dg") -> list[str]:
    return gather_cmd(layout.results_dir, layout.dg_tsv, report=report)
