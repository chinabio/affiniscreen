"""Login-node job throttler for the Amber MD / OpenFE workflow (v2.5.0+).

Standalone, dependency-free (stdlib only). Launched DETACHED via nohup by the
Setup & Launch page ONLY when the user enables the optional submission throttle.
It is never imported or run in the default launch path, so default behaviour is
unchanged.

It reads a JSON job-list (a list of {"cmd": [...], "cwd": "...", "log": "..."} )
and runs them so that no more than --max-inflight of THIS batch's jobs are
"in the queue" (PEND+RUN) at any time. Each job in the list is itself a command
that submits exactly one cluster job (e.g. `bsub < script` or a tool that calls
bsub once) OR runs locally; the throttler tags every submission with a unique
LSF job-name prefix and counts matching names via `bjobs` to decide when to
release the next one.

Because it runs under nohup on the login node, it survives the Streamlit page
being closed -- which the previous GUI-side throttle did not.

Usage:
    python -m amber_md.throttle_submit --joblist jobs.json \
        --max-inflight 8 --name-prefix mmgbsa_run42 [--poll 30] [--submit-gap 3]
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


def _log(msg: str, logf):
    line = f"{datetime.now().isoformat(timespec='seconds')}  {msg}\n"
    logf.write(line)
    logf.flush()


def _count_inflight(name_prefix: str) -> int:
    """Number of this batch's jobs currently PEND or RUN, counted by bjobs.

    Matches on JOB_NAME starting with name_prefix. Returns a large sentinel on
    bjobs failure so the throttler waits rather than flooding the queue.
    """
    user = os.environ.get("USER", "")
    try:
        cp = subprocess.run(
            ["bjobs", "-u", user, "-o", "job_name stat", "-noheader"],
            capture_output=True, text=True, timeout=30)
    except Exception:
        return -1  # signal "unknown"
    if cp.returncode != 0:
        return -1
    n = 0
    for ln in cp.stdout.splitlines():
        parts = ln.split()
        if not parts:
            continue
        jn = parts[0]
        stat = parts[1] if len(parts) > 1 else ""
        if jn.startswith(name_prefix) and stat in ("PEND", "RUN", "PROV", "WAIT"):
            n += 1
    return n


def run(joblist_path, max_inflight, name_prefix, poll, submit_gap, logf):
    jobs = json.loads(Path(joblist_path).read_text())
    total = len(jobs)
    _log(f"throttler start: {total} jobs, max_inflight={max_inflight}, "
         f"prefix={name_prefix}", logf)
    submitted = 0
    for i, job in enumerate(jobs):
        # wait until there is room
        while True:
            inflight = _count_inflight(name_prefix)
            if inflight < 0:
                _log("bjobs unavailable; waiting before next submission", logf)
                time.sleep(poll)
                continue
            if inflight < max_inflight:
                break
            _log(f"inflight={inflight} >= {max_inflight}; waiting {poll}s "
                 f"(submitted {submitted}/{total})", logf)
            time.sleep(poll)

        cmd = job["cmd"]
        cwd = job.get("cwd") or None
        jlog = job.get("log")
        _log(f"submitting job {i+1}/{total}: {' '.join(map(str, cmd))}", logf)
        try:
            if jlog:
                Path(jlog).parent.mkdir(parents=True, exist_ok=True)
                with open(jlog, "a") as jf:
                    jf.write(f"\n# Throttled-submit at {datetime.now().isoformat()}\n")
                    jf.write(f"# Command: {' '.join(map(str, cmd))}\n\n")
                    jf.flush()
                    subprocess.run(cmd, cwd=cwd, stdout=jf,
                                   stderr=subprocess.STDOUT, check=False)
            else:
                subprocess.run(cmd, cwd=cwd, check=False)
            submitted += 1
        except Exception as e:  # noqa: BLE001
            _log(f"ERROR submitting job {i+1}: {e}", logf)
        time.sleep(submit_gap)  # small gap so bjobs reflects the new job

    _log(f"throttler done: submitted {submitted}/{total} jobs", logf)


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m amber_md.throttle_submit")
    p.add_argument("--joblist", required=True,
                   help="JSON file: list of {cmd:[...], cwd, log}")
    p.add_argument("--max-inflight", type=int, default=8, dest="max_inflight")
    p.add_argument("--name-prefix", required=True, dest="name_prefix",
                   help="LSF job-name prefix used to count this batch via bjobs")
    p.add_argument("--poll", type=int, default=30,
                   help="Seconds between queue checks (default 30)")
    p.add_argument("--submit-gap", type=int, default=3, dest="submit_gap",
                   help="Seconds to wait after each submission (default 3)")
    p.add_argument("--logfile", default=None,
                   help="Where to write throttler progress (default: stderr)")
    a = p.parse_args(argv)
    if a.logfile:
        Path(a.logfile).parent.mkdir(parents=True, exist_ok=True)
        logf = open(a.logfile, "a")
    else:
        logf = sys.stderr
    try:
        run(a.joblist, a.max_inflight, a.name_prefix, a.poll, a.submit_gap, logf)
    finally:
        if logf is not sys.stderr:
            logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
