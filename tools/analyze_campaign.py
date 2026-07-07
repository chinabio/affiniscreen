#!/usr/bin/env python3
"""Campaign-wide ABFE analyzer: analyze every leg of an edge and combine to dG_bind.

Re-analyzes all legs under a fep/ edge directory (no LSF, no re-simulation) and
assembles the binding free energy using the SAME thermodynamic cycle and gating
the in-pipeline cycle-closer uses:

    dG_complex = complex_decharge + complex_vdw + complex_restraint + Boresch
    dG_solvent = solvent_decharge + solvent_vdw
    dG_bind    = -(dG_complex - dG_solvent) + dG_charge_correction

Each leg is scored with the same logic as tools/analyze_leg.py:
    * Option-A restraint leg (no clambda) -> analytic Boresch value
    * TI legs (decharge / vdw)            -> FEPAnalyzer / MBAR

Usage
-----
    # from the workflow root (so `amber_md` is importable):
    python tools/analyze_campaign.py <fep_dir>
    python tools/analyze_campaign.py <fep_dir> --temp 298.0
    python tools/analyze_campaign.py <campaign_root> --recurse   # many edges
    python tools/analyze_campaign.py <fep_dir> --json-only

<fep_dir> is the directory that directly contains the leg sub-directories
(complex_decharge/, complex_vdw/, complex_restraint/, solvent_decharge/,
solvent_vdw/).  With --recurse, any directory beneath <campaign_root> that
contains a complex_vdw/ leg is treated as an edge and analyzed.

Writes per-edge ABFE_RESULT.json + ABFE_RESULT.txt into each edge dir.

Exit codes:
    0  every analyzed edge is trusted (complete + reliable + sane)
    1  at least one edge is UNTRUSTED (incomplete / unreliable / unphysical)
    2  no dG_bind could be produced for any edge (missing/failed legs)
    4  usage / discovery error
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# canonical leg names
COMPLEX_LEGS = ["complex_decharge", "complex_vdw", "complex_restraint"]
SOLVENT_LEGS = ["solvent_decharge", "solvent_vdw"]
ALL_LEGS = COMPLEX_LEGS + SOLVENT_LEGS
DG_BIND_SANITY_KCAL = 25.0


def _lambdas_for_leg(leg_dir: Path):
    """Reuse analyze_leg's discovery: analyze_*.lsf first, then lambda_* dirs."""
    for lsf in sorted(leg_dir.glob("analyze_*.lsf")):
        txt = lsf.read_text(errors="ignore")
        m = re.search(r'"([0-9.,\s]+)"\.split\(","\)', txt)
        if m:
            try:
                vals = [float(x) for x in m.group(1).split(",") if x.strip()]
                if vals:
                    return vals
            except ValueError:
                pass
    vals = []
    for d in sorted(leg_dir.glob("lambda_*")):
        m = re.match(r"lambda_([0-9]*\.?[0-9]+)$", d.name)
        if m and d.is_dir():
            vals.append(float(m.group(1)))
    return sorted(set(vals))


def analyze_one_leg(leg_dir: Path, temp: float, force: bool = False):
    """Return the summary dict for a single leg, or None if the leg is absent.

    Resume-safe (mirrors the package's own _run cache): if the leg already has a
    summary.json with a non-null dG_kcal_mol, reuse it instead of re-analyzing,
    unless force=True. This makes finished campaigns fast to combine and avoids
    needlessly re-running MBAR (which needs pandas/alchemlyb) on legs that are
    already scored.
    """
    if not leg_dir.is_dir():
        return None
    sj = leg_dir / "summary.json"
    if not force and sj.exists():
        try:
            cached = json.loads(sj.read_text())
            if cached.get("dG_kcal_mol") is not None:
                return cached
        except Exception:
            pass
    from amber_md.fep import analyze_restraint_leg_optionA, FEPAnalyzer
    lambdas = _lambdas_for_leg(leg_dir)
    if not lambdas:
        return None
    res = analyze_restraint_leg_optionA(leg_dir, lambdas, temperature_K=temp)
    if res is None:
        res = FEPAnalyzer(leg_dir, lambdas, temperature_K=temp).run()
    try:
        sj.write_text(json.dumps(res, indent=2, default=str))
    except Exception:
        pass
    return res


def _g(r):
    return None if r is None else r.get("dG_kcal_mol")


def _leg_ok(r):
    return (r is not None and r.get("dG_kcal_mol") is not None
            and r.get("complete", True))


def analyze_edge(fep_dir: Path, temp: float, json_only: bool, force: bool = False):
    """Analyze all legs of one edge and write/return the combined ABFE result."""
    res = {name: analyze_one_leg(fep_dir / name, temp, force=force) for name in ALL_LEGS}

    cdc, cvd, crest = res["complex_decharge"], res["complex_vdw"], res["complex_restraint"]
    sdc, svd = res["solvent_decharge"], res["solvent_vdw"]

    # charge correction (optional)
    dG_charge = 0.0
    ccf = fep_dir / "charge_correction.json"
    if ccf.exists():
        try:
            dG_charge = float(json.loads(ccf.read_text()).get("dG_charge_correction_kcal_mol", 0.0))
        except Exception:
            dG_charge = 0.0

    # complex leg total = decharge + vdw + restraint contribution.
    # IMPORTANT: avoid double-counting the Boresch term, and support BOTH layouts:
    #
    #   (1) WITH a restraint MD leg (default pre-2.5.68 / opt-in --restraint-leg):
    #       * Option-A restraint leg: analyze_restraint_leg_optionA already sets
    #         dG_kcal_mol = Boresch (dG_uncorrected = 0). Add it ONCE.
    #       * TI restraint leg: dG_kcal_mol is the MD restraint FE and
    #         dG_boresch_correction is a SEPARATE analytic term -> add both.
    #
    #   (2) NO restraint leg (v2.5.68 default --no-restraint-leg): the analytic
    #       Boresch term is written to complex_vdw/boresch_correction.txt instead.
    #       There is no complex_restraint dir; read the term from complex_vdw.
    has_restraint_leg = crest is not None and _g(crest) is not None
    # analytic Boresch term carried on complex_vdw when there is no restraint leg
    vdw_bcorr = None
    vcf = fep_dir / "complex_vdw" / "boresch_correction.txt"
    if not has_restraint_leg and vcf.exists():
        try:
            vdw_bcorr = float(vcf.read_text().strip())
        except Exception:
            vdw_bcorr = None
    ctot = None
    if _g(cdc) is not None and _g(cvd) is not None and has_restraint_leg:
        # layout (1)
        ctot = _g(cdc) + _g(cvd) + _g(crest)
        is_optionA = (crest.get("method") == "analytic_boresch"
                      or crest.get("leg_type") == "restraint_optionA")
        bcorr = crest.get("dG_boresch_correction")
        if bcorr is not None and not is_optionA:
            ctot += bcorr
    elif _g(cdc) is not None and _g(cvd) is not None and vdw_bcorr is not None:
        # layout (2): no restraint leg, analytic term folded onto complex_vdw
        ctot = _g(cdc) + _g(cvd) + vdw_bcorr
    if _g(sdc) is not None and _g(svd) is not None:
        stot = _g(sdc) + _g(svd)
        stot = _g(sdc) + _g(svd)

    dG_bind = None if (ctot is None or stot is None) else -(ctot - stot) + dG_charge
    # gating, mirroring the in-pipeline cycle-closer.
    # When the analytic Boresch term is folded onto complex_vdw (no restraint
    # leg), complex_restraint is legitimately absent -> don't treat it as a
    # missing/incomplete leg.
    required_legs = list(ALL_LEGS)
    if not has_restraint_leg and vdw_bcorr is not None:
        required_legs = [k for k in ALL_LEGS if k != "complex_restraint"]
    incomplete = sorted(k for k in required_legs if not _leg_ok(res[k]))
    unreliable = sorted(k for k in required_legs
                        if res[k] is not None and res[k].get("dG_reliable") is False)
    mbar_failed = sorted(k for k in required_legs
                         if res[k] is not None and res[k].get("only_ti_survived"))
    sane_mag = (dG_bind is not None) and (abs(dG_bind) <= DG_BIND_SANITY_KCAL)
    trusted = (dG_bind is not None and not incomplete and not unreliable
               and not mbar_failed and sane_mag)

    def Fmt(x):
        return "(failed)" if x is None else ("%+8.3f kcal/mol" % x)

    def Cmpl(r):
        if r is None:
            return "MISSING"
        if r.get("dG_kcal_mol") is None:
            return "FAILED"
        if not r.get("complete", True):
            return "INCOMPLETE %s/%s" % (r.get("n_windows"), r.get("n_requested"))
        return "complete"

    lines = ["============== ABFE RESULT (campaign re-analyze) =============="]
    lines += ["  edge: %s" % fep_dir]
    if has_restraint_leg:
        lines += ["  COMPLEX:  decharge=" + Fmt(_g(cdc)) + "  vdw=" + Fmt(_g(cvd))
                  + "  restraint=" + Fmt(_g(crest)) + "  total(+Boresch)=" + Fmt(ctot)]
    else:
        lines += ["  COMPLEX:  decharge=" + Fmt(_g(cdc)) + "  vdw=" + Fmt(_g(cvd))
                  + "  Boresch(analytic,on vdw)=" + Fmt(vdw_bcorr)
                  + "  total=" + Fmt(ctot)]
    lines += ["  SOLVENT:  decharge=" + Fmt(_g(sdc)) + "  vdw=" + Fmt(_g(svd))
              + "  total=" + Fmt(stot)]
    if dG_charge:
        lines += ["  charge correction: " + Fmt(dG_charge)]
    lines += ["  leg status:"]
    for name in required_legs:
        lines += ["      %-18s %s" % (name, Cmpl(res[name]))]
    if not has_restraint_leg and vdw_bcorr is not None:
        lines += ["      %-18s %s" % ("complex_restraint",
                  "analytic on complex_vdw (no MD leg)")]
    lines += ["---------------------------------------------------------------"]
    if trusted:
        lines += ["  dG_bind = " + Fmt(dG_bind)]
    else:
        why = []
        if dG_bind is None:
            why.append("missing/failed legs")
        if incomplete:
            why.append("incomplete: " + ", ".join(incomplete))
        if unreliable:
            why.append("reliability gate failed: " + ", ".join(unreliable))
        if mbar_failed:
            why.append("MBAR/BAR failed (TI-only): " + ", ".join(mbar_failed))
        if dG_bind is not None and abs(dG_bind) > DG_BIND_SANITY_KCAL:
            why.append("dG_bind magnitude %.1f > %.0f (unphysical)"
                       % (abs(dG_bind), DG_BIND_SANITY_KCAL))
        lines += ["  dG_bind = " + Fmt(dG_bind) + "   *** UNTRUSTED ("
                  + ("; ".join(why) or "untrusted") + ") ***"]
    lines += ["==============================================================="]
    txt = "\n".join(lines) + "\n"

    out_json = {
        "edge_dir": str(fep_dir),
        "dG_bind_kcal_mol": dG_bind,
        "trusted": trusted,
        "dG_complex_plus_restr_kcal_mol": ctot,
        "dG_solvent_kcal_mol": stot,
        "dG_charge_correction_kcal_mol": dG_charge,
        "incomplete_legs": incomplete,
        "unreliable_legs": unreliable,
        "mbar_failed_legs": mbar_failed,
        "legs": {name: _g(res[name]) for name in ALL_LEGS},
        "complex_restraint_boresch": (None if crest is None
                                      else crest.get("dG_boresch_correction")),
    }
    try:
        (fep_dir / "ABFE_RESULT.txt").write_text(txt)
        (fep_dir / "ABFE_RESULT.json").write_text(json.dumps(out_json, indent=2, default=str))
    except Exception as e:
        print("[analyze_campaign] WARN: could not write ABFE_RESULT.* (%s)" % e,
              file=sys.stderr)

    if not json_only:
        print(txt)
    return out_json


def find_edges(root: Path):
    """Any directory containing a complex_vdw/ sub-dir is treated as an edge."""
    edges = []
    if (root / "complex_vdw").is_dir():
        edges.append(root)
    for d in sorted(root.rglob("complex_vdw")):
        if d.is_dir():
            edges.append(d.parent)
    # dedupe, preserve order
    seen, uniq = set(), []
    for e in edges:
        if e not in seen:
            seen.add(e); uniq.append(e)
    return uniq


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="analyze_campaign.py",
        description="Re-analyze every leg of an ABFE edge (or many edges) and "
                    "combine to dG_bind.")
    ap.add_argument("fep_dir", help="edge fep/ dir (contains the leg sub-dirs), "
                                    "or a campaign root with --recurse")
    ap.add_argument("--temp", "--temperature", dest="temp", type=float, default=298.0)
    ap.add_argument("--recurse", action="store_true",
                    help="treat fep_dir as a campaign root and analyze every edge "
                         "(any dir containing complex_vdw/) beneath it")
    ap.add_argument("--force", action="store_true",
                    help="re-analyze every leg even if a valid summary.json exists")
    ap.add_argument("--json-only", action="store_true",
                    help="suppress per-edge text; still writes ABFE_RESULT.* files")
    ap.add_argument("--csv", default=None,
                    help="write a one-row-per-edge CSV here. With --recurse and no "
                         "--csv, a campaign_summary.csv is written at the root.")
    a = ap.parse_args(argv)
    a = ap.parse_args(argv)

    root = Path(a.fep_dir).expanduser().resolve()
    if not root.is_dir():
        print(f"[analyze_campaign] ERROR: not a directory: {root}", file=sys.stderr)
        return 4

    try:
        import amber_md.fep  # noqa: F401
    except Exception as e:
        print(f"[analyze_campaign] ERROR: cannot import amber_md.fep ({e}).\n"
              f"  Run from the workflow root or set PYTHONPATH.", file=sys.stderr)
        return 4

    edges = find_edges(root) if a.recurse else (
        [root] if (root / "complex_vdw").is_dir() else find_edges(root))
    if not edges:
        print(f"[analyze_campaign] ERROR: no edges found under {root} "
              f"(looked for a complex_vdw/ leg).", file=sys.stderr)
        return 4

    print(f"[analyze_campaign] {len(edges)} edge(s) to analyze; T={a.temp} K")
    results = [analyze_edge(e, a.temp, a.json_only, force=a.force) for e in edges]

    # campaign summary table
    print("\n================= CAMPAIGN SUMMARY =================")
    print("%-50s %14s  %s" % ("edge", "dG_bind", "status"))
    print("-" * 75)
    any_bind = False
    all_trusted = True
    for r in results:
        b = r["dG_bind_kcal_mol"]
        any_bind = any_bind or (b is not None)
        all_trusted = all_trusted and bool(r["trusted"])
        label = Path(r["edge_dir"]).name or r["edge_dir"]
        bind = "(failed)" if b is None else ("%+8.3f" % b)
        status = "TRUSTED" if r["trusted"] else "UNTRUSTED"
        print("%-50s %14s  %s" % (label[:50], bind, status))
    print("===================================================")

    # one-row-per-edge CSV (sortable ranking artifact)
    csv_path = None
    if a.csv:
        csv_path = Path(a.csv).expanduser()
    elif a.recurse:
        csv_path = root / "campaign_summary.csv"
    if csv_path is not None:
        import csv as _csv
        cols = ["edge", "edge_dir", "dG_bind_kcal_mol", "trusted",
                "dG_complex_plus_restr_kcal_mol", "dG_solvent_kcal_mol",
                "dG_charge_correction_kcal_mol",
                "complex_decharge", "complex_vdw", "complex_restraint",
                "solvent_decharge", "solvent_vdw",
                "complex_restraint_boresch",
                "incomplete_legs", "unreliable_legs", "mbar_failed_legs"]
        try:
            csv_path.parent.mkdir(parents=True, exist_ok=True)
            with open(csv_path, "w", newline="") as fh:
                w = _csv.writer(fh)
                w.writerow(cols)
                # rank trusted edges by dG_bind (most negative = tightest binder) first
                def _key(r):
                    b = r.get("dG_bind_kcal_mol")
                    return (0 if r.get("trusted") else 1,
                            b if b is not None else float("inf"))
                for r in sorted(results, key=_key):
                    legs = r.get("legs", {})
                    w.writerow([
                        Path(r["edge_dir"]).name or r["edge_dir"],
                        r["edge_dir"],
                        r.get("dG_bind_kcal_mol"),
                        r.get("trusted"),
                        r.get("dG_complex_plus_restr_kcal_mol"),
                        r.get("dG_solvent_kcal_mol"),
                        r.get("dG_charge_correction_kcal_mol"),
                        legs.get("complex_decharge"), legs.get("complex_vdw"),
                        legs.get("complex_restraint"),
                        legs.get("solvent_decharge"), legs.get("solvent_vdw"),
                        r.get("complex_restraint_boresch"),
                        "|".join(r.get("incomplete_legs", [])),
                        "|".join(r.get("unreliable_legs", [])),
                        "|".join(r.get("mbar_failed_legs", [])),
                    ])
            print(f"[analyze_campaign] CSV written: {csv_path}")
        except Exception as e:
            print(f"[analyze_campaign] WARN: could not write CSV ({e})", file=sys.stderr)

    if not any_bind:
        return 2
    return 0 if all_trusted else 1

if __name__ == "__main__":
    sys.exit(main())