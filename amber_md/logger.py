
# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

import logging
from pathlib import Path
def get_logger(name="amber_md", log_file=None, level=logging.INFO):
    lg = logging.getLogger(name)
    if lg.handlers: return lg
    lg.setLevel(level)
    fmt = logging.Formatter("%(asctime)s [%(levelname)-7s] %(name)s :: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(); sh.setFormatter(fmt); lg.addHandler(sh)
    if log_file:
        log_file = Path(log_file); log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file); fh.setFormatter(fmt); lg.addHandler(fh)
    return lg
