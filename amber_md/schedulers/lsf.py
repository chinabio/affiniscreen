"""LSF backend -- wraps the historical `bsub`/`bjobs`/`bkill` behaviour."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations
from typing import Sequence
from .base import Scheduler


class LSFScheduler(Scheduler):
    name = "lsf"

    def submit_command(self, script, *, job_name, queue, walltime, project="",
                       n_slots=1, gpu=0, extra=()):
        s = self.cfg.scheduler
        cmd = [s.submit_cmd, "-J", job_name, "-q", queue, "-W", walltime,
               "-n", str(n_slots)]
        if project:
            cmd += ["-P", project]
        if gpu:
            # the cluster `gpu`: request GPUs via rusage so cores don't inflate GPUs.
            cmd += ["-R", f"rusage[ngpus_physical={gpu}]"]
        sel = self.resource_select(self.cfg.resources.avoid_hosts)
        if sel:
            cmd += ["-R", sel]
        cmd += list(extra)
        cmd += ["<", script]
        return cmd

    def query_command(self, user=""):
        cmd = [self.cfg.scheduler.query_cmd]
        if user:
            cmd += ["-u", user]
        return cmd

    def cancel_command(self, job_id):
        return [self.cfg.scheduler.cancel_cmd, str(job_id)]

    def resource_select(self, avoid_hosts: Sequence[str] = ()):
        hosts = sorted(h for h in (avoid_hosts or ()) if h)
        if not hosts:
            return ""
        clause = " && ".join(f"hname!={h}" for h in hosts)
        return f"select[{clause}]"
