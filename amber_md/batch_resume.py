"""Aggregate analysis for a completed batch screening (v2.4.0).

Usage:
    python -m amber_md.batch_resume <batch_dir>
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import sys, json, subprocess, argparse
from pathlib import Path
from .logger import get_logger
from .config import WorkflowConfig
from .pipeline import AmberPipeline
log = get_logger()

def _bjobs_status(jobid):
    if jobid in ("FAILED", "TIMEOUT", "UNKNOWN", None, ""):
        return jobid or "UNKNOWN"
    try:
        cp = subprocess.run(["bjobs", "-noheader", "-o", "stat", jobid],
                            capture_output=True, text=True, timeout=10)
        return cp.stdout.strip() or "UNKNOWN"
    except Exception:
        return "ERROR"

def _gbsa_score(workdir):
    f = workdir / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
    if not f.exists(): return None
    for line in f.read_text().splitlines():
        s = line.strip()
        if s.startswith("DELTA TOTAL") and len(s.split()) >= 3:
            try: return float(s.split()[-3])
            except Exception: pass
    return None

def _rmsd_summary(workdir):
    f = workdir / "analysis" / "rmsd.dat"
    if not f.exists(): return None, None
    vals = []
    for line in f.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"): continue
        parts = s.split()
        if len(parts) >= 2:
            try: vals.append(float(parts[1]))
            except Exception: pass
    if not vals: return None, None
    return sum(vals) / len(vals), max(vals)

def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m amber_md.batch_resume",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("batch_dir", type=Path,
        help="Path to batch directory containing batch_manifest.json")
    p.add_argument("--analyze-done", action="store_true", default=True,
        help="Auto-run analysis on jobs whose status is DONE.")
    p.add_argument("--no-analyze", dest="analyze_done", action="store_false")
    a = p.parse_args(argv)

    batch_dir = a.batch_dir.resolve()
    manifest_file = batch_dir / "batch_manifest.json"
    if not manifest_file.exists():
        sys.exit(f"ERROR: {manifest_file} not found.")
    manifest = json.loads(manifest_file.read_text())
    log.info("Loaded %d entries from %s", len(manifest), manifest_file.name)

    rows = []
    for m in manifest:
        name = m["name"]; jobid = m["lsf_jobid"]; workdir = Path(m["workdir"])
        status = _bjobs_status(jobid)
        traj = workdir / "jobs" / "prod.nc"
        gbsa_file = workdir / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
        if (a.analyze_done and status == "DONE" and traj.exists() and not gbsa_file.exists()):
            log.info("[%s] DONE -- running analysis (resume)", name)
            try:
                cfg = WorkflowConfig.load(workdir / "workflow_config.json")
                cfg.work_dir = workdir
                pipe = AmberPipeline(cfg)
                pipe.stage_analyze(workdir / "build" / "complex.prmtop")
            except Exception as e:
                log.error("[%s] analysis failed: %s", name, e)
        mean_rmsd, max_rmsd = _rmsd_summary(workdir)
        gbsa = _gbsa_score(workdir)
        rows.append({
            "name": name, "lsf_jobid": jobid, "status": status,
            "workdir": str(workdir),
            "trajectory_exists": traj.exists(),
            "trajectory_MB": round(traj.stat().st_size/1e6, 1) if traj.exists() else None,
            "rmsd_mean_A": round(mean_rmsd, 2) if mean_rmsd is not None else None,
            "rmsd_max_A":  round(max_rmsd, 2)  if max_rmsd  is not None else None,
            "gbsa_dG_kcal_per_mol": round(gbsa, 2) if gbsa is not None else None,
        })

    csv_path = batch_dir / "batch_summary.csv"
    with open(csv_path, "w") as f:
        cols = ["name","lsf_jobid","status","trajectory_exists","trajectory_MB",
                "rmsd_mean_A","rmsd_max_A","gbsa_dG_kcal_per_mol","workdir"]
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) if r[c] is not None else "" for c in cols) + "\n")

    log.info("\n" + "=" * 70)
    log.info("BATCH SUMMARY -> %s", csv_path)
    log.info("=" * 70)
    log.info("%-20s %-10s %-8s %12s %10s %10s",
             "name", "jobid", "status", "rmsd_mean_A", "rmsd_max", "dG_GBSA")
    log.info("-" * 70)
    for r in rows:
        log.info("%-20s %-10s %-8s %12s %10s %10s",
                 r["name"][:20], r["lsf_jobid"], r["status"],
                 f"{r['rmsd_mean_A']:.2f}" if r["rmsd_mean_A"] is not None else "-",
                 f"{r['rmsd_max_A']:.2f}"  if r["rmsd_max_A"]  is not None else "-",
                 f"{r['gbsa_dG_kcal_per_mol']:+.2f}" if r["gbsa_dG_kcal_per_mol"] is not None else "-")
    by_status = {}
    for r in rows:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
    log.info("\nStatus tally: %s",
             "  ".join(f"{k}={v}" for k, v in sorted(by_status.items())))
    return 0

if __name__ == "__main__":
    sys.exit(main() or 0)
