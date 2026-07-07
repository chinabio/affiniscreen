"""SLURM backend -- STUB.

The command shapes below are the intended mapping, but this backend is NOT yet
validated on a real SLURM cluster. Every method raises SchedulerError so a
misconfigured site fails loudly instead of silently generating broken jobs.
Remove the guards and test end-to-end before enabling `scheduler.type: slurm`.
"""
from __future__ import annotations
from typing import Sequence
from .base import Scheduler, SchedulerError

_MSG = ("SLURM support is a stub and not yet validated. Set scheduler.type "
        "to 'lsf' in your site_config, or implement + test amber_md/schedulers/"
        "slurm.py for your cluster.")


class SlurmScheduler(Scheduler):
    name = "slurm"

    def submit_command(self, script, *, job_name, queue, walltime, project="",
                       n_slots=1, gpu=0, extra=()):
        # Intended mapping (for the future implementer):
        #   sbatch -J <job> -p <partition> -t <walltime> -n <slots>
        #          [-A <account>] [--gres=gpu:<n>] <script>
        raise SchedulerError(_MSG)

    def query_command(self, user=""):
        raise SchedulerError(_MSG)

    def cancel_command(self, job_id):
        raise SchedulerError(_MSG)

    def resource_select(self, avoid_hosts: Sequence[str] = ()):
        # SLURM analog: --exclude=<host1,host2>
        raise SchedulerError(_MSG)
