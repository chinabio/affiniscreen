"""
amber_config.py  --  Build a valid WorkflowConfig from wizard parameters
(v2.5.0, Phase 4b).

The Setup & Launch wizard collects advanced Amber settings (igb, MM-GBSA frame
window/stride/decomposition, FEP lambda schedule, masks, protonation) that
run_amber.py only honours via its --config file (NOT via CLI flags). This helper
starts from the canonical dataclass DEFAULTS, overrides only what the wizard
exposes, and uses WorkflowConfig.save() so the written file is guaranteed to
load back through WorkflowConfig.load(). Returns the config path for `--config`.

Boresch restraints are NOT a run_amber.py concept; they belong to the
fep_driver ABFE entrypoint (--auto-boresch / --boresch-json). Those flags are
emitted by the wizard's command builder, not here.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

from pathlib import Path

# Canonical defaults: import the real dataclasses so we never drift from them.
from amber_md.config import (
    WorkflowConfig, MMGBSAConfig, FEPConfig, FEPWorkflowConfig,
)

DEFAULT_LAMBDAS = FEPConfig().lambdas


# ---------------------------------------------------------------------------
# Parsers (pure; raise ValueError with human messages on bad input)
# ---------------------------------------------------------------------------
def parse_lambda_string(text: str) -> list[float] | None:
    if not text or not text.strip():
        return None
    raw = text.replace(",", " ").split()
    try:
        vals = [float(x) for x in raw]
    except ValueError as e:
        raise ValueError(f"Non-numeric lambda value: {e}")
    if len(vals) < 2:
        raise ValueError("Need at least 2 lambda windows (0.0 and 1.0).")
    if any(v < 0.0 or v > 1.0 for v in vals):
        raise ValueError("Lambda values must be within [0, 1].")
    if any(b <= a for a, b in zip(vals, vals[1:])):
        raise ValueError("Lambda values must be strictly increasing.")
    if abs(vals[0]) > 1e-9 or abs(vals[-1] - 1.0) > 1e-9:
        raise ValueError("Lambda schedule must start at 0.0 and end at 1.0.")
    return vals


def parse_protonation_overrides(text: str) -> dict:
    out: dict = {}
    if not text or not text.strip():
        return out
    for tok in text.replace("\n", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" not in tok:
            raise ValueError(f"Bad protonation token '{tok}' (expected KEY:STATE).")
        key, state = (s.strip() for s in tok.split(":", 1))
        if not key or not state:
            raise ValueError(f"Bad protonation token '{tok}'.")
        out[key] = state
    return out


# ---------------------------------------------------------------------------
# Build a WorkflowConfig from wizard params, overriding only exposed fields.
# ---------------------------------------------------------------------------
def build_config(P: dict, *, work_dir: str, method: str,
                 protein: str | None = None, ligand: str | None = None,
                 lig_resname: str = "LIG") -> WorkflowConfig:
    cfg = WorkflowConfig()                      # canonical defaults
    cfg.work_dir = Path(work_dir)
    cfg.ligand_resname = lig_resname
    cfg.ligand_input = Path(ligand) if ligand else None
    cfg.protein_input = Path(protein) if protein else None

    # --- MD simulation time (final52): honor the GUI Production/Equil fields.
    # Previously these were dropped, so a --config launch silently used the
    # MDConfig default (was 50 ns) regardless of what the GUI showed.
    _prod_ns = P.get("prod_ns", P.get("complex_ns", None))
    _equil_ns = P.get("equil_ns", None)
    if _prod_ns is not None:
        cfg.md.prod_nsteps = int(round(float(_prod_ns) * 1e6 / 1))   # dt=1 fs (v2.5.62)
    if _equil_ns is not None:
        cfg.md.equil_nsteps = int(round(float(_equil_ns) * 1e6 / 1))  # dt=1 fs (v2.5.62)

    # --- MM-GBSA ---
    cfg.mmgbsa = MMGBSAConfig(
        enabled=(method == "MM-GBSA"),
        igb=int(P.get("igb", MMGBSAConfig().igb)),
        salt_conc=float(P.get("salt_conc", MMGBSAConfig().salt_conc)),
        start_frame=int(P.get("mmgbsa_start_frame", 1)),
        end_frame=int(P.get("mmgbsa_end_frame", 0)),
        stride=int(P.get("mmgbsa_stride", 1)),
        decomposition=bool(P.get("decomposition", False)),
        decomp_residues=str(P.get("decomp_residues", "")),
    )

    # --- FEP params ---
    lambdas = P.get("lambdas")
    if isinstance(lambdas, str):
        lambdas = parse_lambda_string(lambdas)
    if not lambdas:
        lambdas = list(DEFAULT_LAMBDAS)
    fp = FEPConfig()
    fp.lambdas = tuple(lambdas)
    if method == "ABFE":
        # final55 FIX: FEPConfig defaults to RBFE-style two-ligand masks
        # (:LIA / :LIB), which are WRONG for an absolute (single-topology)
        # decoupling and reference residues that don't exist in an ABFE prmtop.
        # Mirror what fep_driver._derive_abfe_masks() actually does so this
        # provenance config matches the run: region1 = :<resname>, region2 = "".
        _lig = ":" + str(lig_resname).lstrip(":")
        fp.timask1 = P.get("timask1", _lig)
        fp.timask2 = P.get("timask2", "")
        fp.scmask1 = P.get("scmask1", _lig)
        fp.scmask2 = P.get("scmask2", "")
        fp.crgmask = P.get("crgmask", _lig)
    else:
        fp.timask1 = P.get("timask1", fp.timask1)
        fp.timask2 = P.get("timask2", fp.timask2)
        fp.scmask1 = P.get("scmask1", P.get("timask1", fp.scmask1))
        fp.scmask2 = P.get("scmask2", P.get("timask2", fp.scmask2))
    fp.temperature_K = float(P.get("temperature_K", fp.temperature_K))
    fp.cutoff_A = float(P.get("cutoff_A", fp.cutoff_A))
    fp.use_fine_restraint_lambdas = bool(P.get("use_fine_restraint_lambdas", False))
    cfg.fep = FEPWorkflowConfig(enabled=method in ("ABFE", "RBFE"), params=fp)

    # --- protonation ---
    cfg.auto_protonation = bool(P.get("auto_protonation", True))
    overrides = parse_protonation_overrides(P.get("protonation_overrides", ""))
    cfg.protonation_overrides = overrides or None

    cfg.submit = bool(P.get("submit", True))
    cfg.monitor = bool(P.get("monitor", False))
    return cfg


def write_config(P: dict, *, work_dir: str, method: str, **kw) -> Path:
    """Build + save the config JSON under work_dir; return its path."""
    cfg = build_config(P, work_dir=work_dir, method=method, **kw)
    wd = Path(work_dir).expanduser()
    wd.mkdir(parents=True, exist_ok=True)
    out = wd / "wizard_config.json"
    cfg.save(out)                              # canonical serializer
    return out
