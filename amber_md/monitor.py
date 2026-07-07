
# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import time
from .utils import run
from .logger import get_logger
log = get_logger()
TERMINAL = {"DONE","EXIT","ZOMBI","UNKWN"}

class JobMonitor:
    def __init__(self, poll_interval_s=60, max_wait_s=7*24*3600):
        self.poll = poll_interval_s; self.max = max_wait_s
    def status(self, job_id):
        cp = run(["bjobs","-noheader","-o","stat",job_id], check=False, capture=True)
        st = (cp.stdout or "").strip().split()
        return st[0] if st else "UNKWN"
    def wait(self, job_id):
        elapsed=0; last=""
        while elapsed < self.max:
            st = self.status(job_id)
            if st != last: log.info("Job %s status: %s", job_id, st); last = st
            if st in TERMINAL: return st
            time.sleep(self.poll); elapsed += self.poll
        return "TIMEOUT"
