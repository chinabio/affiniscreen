"""v2.5.15 regression: window recovery (purge + clean rebuild)."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import shutil
from pathlib import Path
from amber_md.abfe_self_heal_cli import (
    recovery_stages, purge_poisoned_restarts, _clean_start_coord,
    _PHYS_BLOWUP_CLASSES)


def _detonated(root, with_dens=True):
    wd = Path(root)
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    (wd / "system.inpcrd").write_text("CLEAN\n")
    (wd / "system.prmtop").write_text("P\n")
    (wd / "min.rst").write_text("CLEAN MIN\n")
    if with_dens:
        (wd / "dens.in").write_text("&cntrl\n dt=0.001,\n/\n")
        (wd / "dens.rst").write_text("DETONATED\n")
    (wd / "eq.in").write_text("&cntrl\n dt=0.002,\n/\n")
    (wd / "eq.rst").write_text("DETONATED\n")
    (wd / "prod.in").write_text("&cntrl\n irest=1, dt=0.002,\n/\n")
    (wd / "prod.out").write_text("vlimit exceeded\n")
    return wd


def test_blowup_full_redo_from_clean_coords(tmp_path):
    wd = _detonated(tmp_path / "win")
    seq, full = recovery_stages(wd, "blowup_temperature")
    assert full is True
    assert [s[0] for s in seq] == ["dens", "eq", "prod"]
    assert seq[0][2] == "min.rst"          # clean start
    removed = purge_poisoned_restarts(wd, full)
    assert "eq.rst" in removed and "dens.rst" in removed and "prod.out" in removed
    assert not (wd / "eq.rst").exists()


def test_external_kill_keeps_good_eq_rst(tmp_path):
    wd = Path(tmp_path / "win2"); wd.mkdir()
    (wd / "min.rst").write_text("CLEAN\n")
    (wd / "eq.rst").write_text("GOOD\n")
    (wd / "prod.in").write_text("&cntrl\n irest=1,\n/\n")
    (wd / "prod.out").write_text("Terminated\n")
    seq, full = recovery_stages(wd, "external_kill")
    assert full is False
    assert [s[0] for s in seq] == ["prod"] and seq[0][2] == "eq.rst"
    removed = purge_poisoned_restarts(wd, full)
    assert removed == ["prod.out"] and (wd / "eq.rst").exists()


def test_all_phys_classes_full_redo(tmp_path):
    for i, fc in enumerate(sorted(_PHYS_BLOWUP_CLASSES)):
        wd = _detonated(tmp_path / ("w%d" % i))
        _seq, full = recovery_stages(wd, fc)
        assert full is True, fc
