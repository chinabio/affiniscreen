"""
results_lib.py  --  Streamlit-free result parsers shared by the Results pages
(v2.5.0, Phase 3).

Extracted from the original 8_Results.py so both the single-molecule detail view
and the multi-molecule compare/rank view reuse ONE implementation. Pure
functions only -- no `import streamlit` here, so they are unit-testable and
importable from anywhere.

Covered:
  * MM-GBSA          -> mmgbsa_status() / parse_mmpbsa_dat (delegates to batch_aggregate)
  * Amber FEP legs   -> parse_fep_leg(), fep_run_results(), fep_headline()
  * OpenFE ABFE      -> parse_openfe_result()  (the protocol_result JSON)
  * dir discovery    -> list_ligand_subdirs()
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

_SANE_MAX = 500.0


# ---------------------------------------------------------------------------
# Sanity + headline-row selection (Amber FEP estimator CSVs)
# ---------------------------------------------------------------------------
def is_sane_row(row) -> bool:
    v = row.get("dG_kcal_mol")
    try:
        if v is None or not math.isfinite(float(v)):
            return False
        if abs(float(v)) > _SANE_MAX:
            return False
    except (TypeError, ValueError):
        return False
    err = row.get("err_kcal_mol")
    if err is not None and pd.notna(err):
        try:
            if not math.isfinite(float(err)):
                return False
        except (TypeError, ValueError):
            return False
    return True


def pick_headline_row(df: "pd.DataFrame"):
    df = df.copy()
    df["dG_kcal_mol"] = pd.to_numeric(df["dG_kcal_mol"], errors="coerce")
    df["err_kcal_mol"] = pd.to_numeric(df["err_kcal_mol"], errors="coerce")
    if "headline" in df.columns:
        marked = df[df["headline"].astype(str).str.contains(r"\*", regex=True,
                                                            na=False)]
        if not marked.empty:
            return marked.iloc[0]
    for est in ("TI", "BAR", "MBAR", "TI_legacy"):
        sub = df[df.estimator == est]
        if not sub.empty and is_sane_row(sub.iloc[0]):
            return sub.iloc[0]
    return None


# ---------------------------------------------------------------------------
# Amber FEP legs
# ---------------------------------------------------------------------------
def parse_fep_leg(leg_dir: Path) -> dict | None:
    """Prefer summary.json (analyzer-canonical), fall back to dG_estimators.csv."""
    leg_dir = Path(leg_dir)
    out: dict = {"leg": leg_dir.name}

    sj = leg_dir / "summary.json"
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            out["dG_kcal_mol"] = s.get("dG_kcal_mol")
            out["err_kcal_mol"] = (s.get("estimators", {})
                                   .get(s.get("headline_estimator", ""), {})
                                   .get("err_kcal_mol"))
            out["estimator"] = s.get("headline_estimator")
            out["uncorrected"] = s.get("dG_uncorrected_kcal_mol")
            out["boresch_correction"] = s.get("dG_boresch_correction")
            out["source"] = "summary.json"
        except Exception:
            pass

    est_csv = leg_dir / "dG_estimators.csv"
    if est_csv.exists():
        try:
            df = pd.read_csv(est_csv)
            out["estimators_df"] = df
            if "dG_kcal_mol" not in out or out.get("dG_kcal_mol") is None:
                row = pick_headline_row(df)
                if row is not None:
                    out["dG_kcal_mol"] = float(row["dG_kcal_mol"])
                    out["err_kcal_mol"] = (float(row["err_kcal_mol"])
                                           if pd.notna(row["err_kcal_mol"])
                                           else None)
                    out["estimator"] = row["estimator"]
                    out["source"] = "dG_estimators.csv"
        except Exception:
            pass

    return out if len(out) > 1 else None


def fep_run_results(fep_root: Path) -> dict:
    fep_root = Path(fep_root)
    results = {}
    for leg in ("complex", "solvent", "absolute"):
        leg_dir = fep_root / leg
        if leg_dir.exists():
            r = parse_fep_leg(leg_dir)
            if r:
                results[leg] = r
    return results


def fep_headline(fep_root: Path) -> tuple[float | None, str, dict]:
    """Headline DG/DDG with the correct formula and provenance label."""
    fep_root = Path(fep_root)
    rj = fep_root / "ABFE_RESULT.json"
    if rj.exists():
        try:
            d = json.loads(rj.read_text())
            return (d.get("dG_bind_kcal_mol"),
                    "ABFE DG_bind (cycle-closed)",
                    {"absolute": {"dG_kcal_mol": d.get("dG_complex_plus_restr_kcal_mol"),
                                  "estimator": d.get("absolute_headline")},
                     "solvent": {"dG_kcal_mol": d.get("dG_solvent_kcal_mol"),
                                 "estimator": d.get("solvent_headline")}})
        except Exception:
            pass

    legs = fep_run_results(fep_root)
    has_cmp = "complex" in legs
    has_abs = "absolute" in legs
    has_sol = "solvent" in legs

    if has_cmp and has_sol:
        try:
            ddg = legs["complex"]["dG_kcal_mol"] - legs["solvent"]["dG_kcal_mol"]
            return ddg, "RBFE DDG (relative)", legs
        except Exception:
            pass
    if has_abs and has_sol:
        try:
            dgb = -(legs["absolute"]["dG_kcal_mol"]
                    - legs["solvent"]["dG_kcal_mol"])
            return dgb, "ABFE DG_bind", legs
        except Exception:
            pass
    if has_abs:
        return (legs["absolute"].get("dG_kcal_mol"),
                "absolute leg only", legs)
    if has_cmp:
        return legs["complex"].get("dG_kcal_mol"), "complex leg only", legs
    if has_sol:
        return legs["solvent"].get("dG_kcal_mol"), "solvent leg only", legs
    return None, "no results", {}


# ---------------------------------------------------------------------------
# OpenFE absolute-binding protocol_result JSON
# ---------------------------------------------------------------------------
def parse_openfe_result(result_json: Path) -> dict | None:
    """Parse an OpenFE *_result.json (AbsoluteBindingProtocolResult etc.).

    Returns a dict with estimate/uncertainty (kcal/mol) and, when available,
    per-leg unit estimates + simple convergence diagnostics. Robust to the
    gufe custom-encoded number dicts ({'magnitude':..,'unit':..}).
    """
    result_json = Path(result_json)
    if not result_json.exists():
        return None
    try:
        data = json.loads(result_json.read_text())
    except Exception:
        return None

    def _mag(x):
        if isinstance(x, dict):
            return x.get("magnitude")
        return x

    out = {"estimate_kcal_mol": _mag(data.get("estimate")),
           "uncertainty_kcal_mol": _mag(data.get("uncertainty")),
           "legs": {}}

    ur = data.get("unit_results", {})
    items = ur.values() if isinstance(ur, dict) else (ur or [])
    for u in items:
        if not isinstance(u, dict):
            continue
        name = u.get("name", "")
        if "Analysis" not in name:
            continue
        leg = "complex" if "complex" in name else (
              "solvent" if "solvent" in name else "other")
        outs = u.get("outputs", {})
        out["legs"][leg] = {
            "unit_estimate": _mag(outs.get("unit_estimate")),
            "unit_estimate_error": _mag(outs.get("unit_estimate_error")),
            "standard_state_correction": _mag(
                outs.get("standard_state_correction")),
            "mbar_overlap_scalar": (outs.get("unit_mbar_overlap", {}) or {}
                                    ).get("scalar"),
            "production_iterations": outs.get("production_iterations"),
            "equilibration_iterations": outs.get("equilibration_iterations"),
        }
    return out


# ---------------------------------------------------------------------------
# MM-GBSA
# ---------------------------------------------------------------------------
def parse_mmpbsa_dat(path: Path):
    """Delegate to the canonical parser in batch_aggregate (single source)."""
    try:
        from amber_md.batch_aggregate import parse_mmpbsa_dat as _p
        return _p(Path(path))
    except Exception:
        return None


def mmgbsa_status(wd: Path):
    """Return (status_text, delta_total_or_None, full_result_or_None)."""
    wd = Path(wd)
    dat = wd / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
    prod = wd / "jobs" / "prod.nc"
    result = parse_mmpbsa_dat(dat)
    if result:
        return "DONE", result.get("delta_total"), result
    if prod.exists():
        return "MD done, MMGBSA pending", None, None
    if (wd / "build" / "complex.prmtop").exists():
        return "BUILD done, MD pending", None, None
    return "Not started / failed", None, None



# ---------------------------------------------------------------------------
# Engine detection for an MM-GBSA ligand dir (final63)
# ---------------------------------------------------------------------------
def mmgbsa_engine(wd: Path) -> str:
    """Best-effort engine label for an MM-GBSA ligand dir: 'OpenMM' or 'Amber'.

    Priority:
      1. mmgbsa/engine.json with {"engine": "..."}  (written by the OpenMM
         runner since final63 -- authoritative & self-describing).
      2. OpenMM-exclusive MD artifacts (production.nc / md_log.txt) for runs
         finished BEFORE final63, which the Amber pmemd path never writes.
      3. default 'Amber'.
    """
    wd = Path(wd)
    marker = wd / "mmgbsa" / "engine.json"
    if marker.exists():
        try:
            eng = (json.loads(marker.read_text()) or {}).get("engine")
            if eng:
                return str(eng)
        except Exception:
            pass
    if (wd / "production.nc").exists() or (wd / "md_log.txt").exists():
        return "OpenMM"
    return "Amber"

# ---------------------------------------------------------------------------
# Directory discovery
# ---------------------------------------------------------------------------
def list_ligand_subdirs(parent: Path) -> list[Path]:
    parent = Path(parent)
    BLOCKLIST = {"logs", "build", "build_solvent", "ligand", "prep",
                 "jobs", "mmgbsa", "fep"}
    by_pattern = list(parent.glob("lig_*")) + list(parent.glob("fep_lig_*"))
    by_content = []
    if parent.is_dir():
        for child in parent.iterdir():
            if not child.is_dir() or child.name in BLOCKLIST:
                continue
            if child.name.startswith("."):
                continue
            if (child / "fep").is_dir() or (child / "mmgbsa").is_dir():
                by_content.append(child)
    seen, ordered = set(), []
    for d in by_pattern + by_content:
        if d.is_dir() and d.name not in seen:
            seen.add(d.name)
            ordered.append(d)
    return sorted(ordered, key=lambda p: p.name)
