"""Scheduler abstraction.

Today the workflow generates and submits LSF (`bsub`) jobs directly in several
places. This package introduces a small, documented interface so a future
SLURM backend can be added in ONE place rather than scattered `if lsf/slurm`
checks. The LSF backend wraps the existing behaviour; the SLURM backend is a
deliberate stub that raises until validated on a real cluster.
"""
from __future__ import annotations
from .base import Scheduler, SchedulerError
from .lsf import LSFScheduler
from .slurm import SlurmScheduler


def get_scheduler(site_cfg=None) -> Scheduler:
    """Return a Scheduler instance for the configured backend."""
    from amber_md import site_config as sc
    cfg = site_cfg or sc.get()
    kind = (cfg.scheduler.type or "lsf").lower()
    if kind == "lsf":
        return LSFScheduler(cfg)
    if kind == "slurm":
        return SlurmScheduler(cfg)
    if kind == "local":
        # Local execution reuses the LSF wrapper's spawn path minus bsub;
        # kept minimal here -- callers that need it can special-case "local".
        return LSFScheduler(cfg)
    raise SchedulerError(f"Unknown scheduler type: {cfg.scheduler.type!r}")
