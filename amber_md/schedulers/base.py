"""Scheduler interface."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Sequence


class SchedulerError(RuntimeError):
    pass


class Scheduler(ABC):
    """Minimal interface every backend must provide.

    The goal is to isolate scheduler-specific syntax (submit/query/cancel and
    the resource-request string) so the rest of the codebase is portable.
    """

    name: str = "base"

    def __init__(self, site_cfg):
        self.cfg = site_cfg

    @abstractmethod
    def submit_command(self, script: str, *, job_name: str, queue: str,
                       walltime: str, project: str = "", n_slots: int = 1,
                       gpu: int = 0, extra: Sequence[str] = ()) -> list[str]:
        """Return the argv used to SUBMIT `script` (e.g. ['bsub', ...])."""

    @abstractmethod
    def query_command(self, user: str = "") -> list[str]:
        """Return the argv used to LIST the user's jobs."""

    @abstractmethod
    def cancel_command(self, job_id: str) -> list[str]:
        """Return the argv used to CANCEL a job."""

    @abstractmethod
    def resource_select(self, avoid_hosts: Sequence[str] = ()) -> str:
        """Return the scheduler-specific host-exclusion / select directive
        (may be empty)."""
