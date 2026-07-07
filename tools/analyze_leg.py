#!/usr/bin/env python3
"""Stand-alone analyzer for a single ABFE/RBFE leg directory (v2.5.65).

Re-analyze a leg that has already finished running -- no LSF, no re-simulation.
Works for BOTH:
  * Option-A restraint legs (plain NPT MD, no clambda) -> analytic Boresch value
  * TI legs (decharge / vdw / solvent)                 -> FEPAnalyzer / MBAR

It runs the SAME logic the pipeline's generated analyze_<leg>.lsf uses, so the
summary.json it writes is identical to what a successful in-pipeline analyze job
would have produced.

Usage
-----
    # from the workflow root (so `amber_md` is importable):
    python tools/analyze_leg.py <leg_dir>
    python tools/analyze_leg.py <leg_dir> --lambdas 0.0,0.15,0.3,1.0
    python tools/analyze_leg.py <leg_dir> --temp 298.0 --json-only

Lambda discovery (when --lambdas is omitted), in order of preference:
    1. the lam_csv embedded in <leg_dir>/analyze_*.lsf  (exactly what was run)
    2. the lambda_<x.xxx> sub-directories present on disk

Exit codes (match the in-pipeline analyzer):
    0  OK (complete, dG produced)
    1  INCOMPLETE (some windows missing; result unreliable)
    2  FAILED (no dG produced)
    3  analyzer raised an exception
    4  usage / discovery error (bad leg dir, no lambdas found)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def _lambdas_from_lsf(leg_dir: Path):
    """Recover the lambda schedule from a generated analyze_*.lsf, if present.

    The analyze LSF embeds:  lambdas = [float(x) for x in "0.000000,0.150000,...".split(",")]
    """
    for lsf in sorted(leg_dir.glob("analyze_*.lsf")):
        txt = lsf.read_text(errors="ignore")
        m = re.search(r'split\(","\)\s*\]\s*', txt)  # locate the lam_csv line
        m = re.search(r'"([0-9.,\s]+)"\.split\(","\)', txt)
        if m:
            try:
                vals = [float(x) for x in m.group(1).split(",") if x.strip()]
                if vals:
                    return vals, f"analyze lsf ({lsf.name})"
            except ValueError:
                pass
    return None, None


def _lambdas_from_dirs(leg_dir: Path):
    """Fall back to scanning lambda_<x.xxx> sub-directories."""
    vals = []
    for d in sorted(leg_dir.glob("lambda_*")):
        if not d.is_dir():
            continue
        m = re.match(r"lambda_([0-9]*\.?[0-9]+)$", d.name)
        if m:
            vals.append(float(m.group(1)))
    vals = sorted(set(vals))
    return (vals, "lambda_* directories") if vals else (None, None)


def discover_lambdas(leg_dir: Path):
    vals, src = _lambdas_from_lsf(leg_dir)
    if vals:
        return vals, src
    return _lambdas_from_dirs(leg_dir)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="analyze_leg.py",
        description="Stand-alone re-analyzer for one ABFE/RBFE leg directory.",
    )
    ap.add_argument("leg_dir", help="leg directory (contains lambda_*/prod.out)")
    ap.add_argument("--lambdas", default=None,
                    help="comma-separated lambda list; if omitted, auto-discovered "
                         "from analyze_*.lsf then from lambda_* dirs")
    ap.add_argument("--temp", "--temperature", dest="temp", type=float, default=298.0,
                    help="temperature in K (default 298.0)")
    ap.add_argument("--json-only", action="store_true",
                    help="suppress the human-readable summary; still writes summary.json")
    a = ap.parse_args(argv)

    leg_dir = Path(a.leg_dir).expanduser().resolve()
    if not leg_dir.is_dir():
        print(f"[analyze_leg] ERROR: not a directory: {leg_dir}", file=sys.stderr)
        return 4

    # import after arg parsing so --help works without the package on path
    try:
        from amber_md.fep import analyze_restraint_leg_optionA, FEPAnalyzer
    except Exception as e:  # pragma: no cover
        print(f"[analyze_leg] ERROR: cannot import amber_md.fep ({e}).\n"
              f"  Run from the workflow root or set PYTHONPATH to it.", file=sys.stderr)
        return 4

    if a.lambdas:
        try:
            lambdas = [float(x) for x in a.lambdas.split(",") if x.strip()]
        except ValueError:
            print(f"[analyze_leg] ERROR: bad --lambdas: {a.lambdas!r}", file=sys.stderr)
            return 4
        src = "--lambdas"
    else:
        lambdas, src = discover_lambdas(leg_dir)
        if not lambdas:
            print(f"[analyze_leg] ERROR: could not discover lambdas in {leg_dir} "
                  f"(no analyze_*.lsf and no lambda_* dirs). Pass --lambdas.",
                  file=sys.stderr)
            return 4

    if not a.json_only:
        print(f"[analyze_leg] leg     : {leg_dir}")
        print(f"[analyze_leg] lambdas : {len(lambdas)} (from {src})")
        print(f"[analyze_leg] temp    : {a.temp} K")

    # Same logic as the generated analyze_<leg>.lsf: Option-A first, MBAR fallback.
    try:
        res = analyze_restraint_leg_optionA(leg_dir, lambdas, temperature_K=a.temp)
        if res is None:
            res = FEPAnalyzer(leg_dir, lambdas, temperature_K=a.temp).run()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[analyze_leg] ANALYZER FAILED: {e}", file=sys.stderr)
        return 3

    (leg_dir / "summary.json").write_text(json.dumps(res, indent=2, default=str))

    if not a.json_only:
        print(f"======== {leg_dir.name} ========")
        print("method / estimator :", res.get("method") or res.get("estimator_used"))
        print("dG (kcal/mol)      :", res.get("dG_kcal_mol"))
        print("windows            : %s/%s" % (res.get("n_windows"), res.get("n_requested")))
        if res.get("missing_windows"):
            print("missing windows    :", res.get("missing_windows"))
        if res.get("dG_boresch_correction") is not None:
            print("  uncorrected dG   :", res.get("dG_uncorrected_kcal_mol"))
            print("  Boresch corr.    :", res.get("dG_boresch_correction"))
        print("summary.json       :", leg_dir / "summary.json")

    if res.get("dG_kcal_mol") is None:
        if not a.json_only:
            print("ANALYZER STATUS: FAILED (no dG produced)")
        return 2
    if not res.get("complete", True):
        if not a.json_only:
            print("ANALYZER STATUS: INCOMPLETE (result UNRELIABLE)")
        return 1
    if not a.json_only:
        print("ANALYZER STATUS: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
