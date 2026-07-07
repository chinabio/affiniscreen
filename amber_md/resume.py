"""Resume analysis on an already-completed MD job.

Usage:
    python -m amber_md.resume <workdir>

This:
  1. Loads {workdir}/workflow_config.json saved by the original run
  2. Skips PREP / BUILD / MD stages
  3. Runs only Stage 4 (analysis) and Stage 5 (FEP if configured)

Useful when you submitted with --no-monitor and want to process the
trajectory later, without rerunning the whole pipeline.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import sys, argparse
from pathlib import Path
from .config import WorkflowConfig
from .pipeline import AmberPipeline
from .logger import get_logger

def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m amber_md.resume",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("workdir", type=Path,
                   help="Path to the original --workdir (must contain workflow_config.json)")
    p.add_argument("--no-gbsa", action="store_true",
                   help="Skip MM-GBSA even if originally enabled.")
    p.add_argument("--no-fep", action="store_true",
                   help="Skip Stage 5 FEP even if originally enabled.")
    a = p.parse_args(argv)

    workdir = Path(a.workdir).resolve()
    cfg_file = workdir / "workflow_config.json"
    if not cfg_file.exists():
        sys.exit(f"ERROR: {cfg_file} not found. Did the original run reach Stage 1?")

    cfg = WorkflowConfig.load(cfg_file)
    # Override work_dir in case the config was loaded from elsewhere
    cfg.work_dir = workdir
    if a.no_gbsa:
        cfg.mmgbsa.enabled = False
    if a.no_fep:
        cfg.fep.enabled = False

    log = get_logger()
    log.info("=" * 64)
    log.info("RESUMING analysis-only run in %s", workdir)
    log.info("=" * 64)

    # Locate the prmtop
    prmtop = workdir / "build" / "complex.prmtop"
    if not prmtop.exists():
        sys.exit(f"ERROR: {prmtop} not found. Original BUILD stage may have failed.")
    log.info("Using prmtop: %s", prmtop)

    # Locate trajectory
    traj = workdir / "jobs" / "prod.nc"
    if not traj.exists():
        sys.exit(
            f"ERROR: {traj} not found.\n"
            f"  Has the GPU job finished yet? Check with: bjobs -u $USER\n"
            f"  Or look in {workdir}/jobs/ for *.out files to see what happened.")
    log.info("Using trajectory: %s (%.1f MB)", traj, traj.stat().st_size / 1e6)

    # Run only the analysis + FEP stages
    p = AmberPipeline(cfg)
    p.stage_analyze(prmtop)
    p.stage_fep()
    log.info("Resume complete.")

if __name__ == "__main__":
    main()
