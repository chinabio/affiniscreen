"""amber_md/abfe_self_heal_cli.py -- run the self-heal loop over a finished
(or partially finished) ABFE run directory.

For each of the four legs it walks every lambda window, and for any window whose
prod.out is missing/incomplete it diagnoses the cause, edits that window's mdin
parameters, and RE-RUNS the window until prod.out is complete (or attempts are
exhausted). Re-run is done by resubmitting just that one window to LSF via the
same single-window mechanism the resume path uses; pass --local to run
pmemd.cuda directly instead (e.g. on an interactive GPU node), or --dry-run to
only diagnose + edit without running.

Usage
-----
    python -m amber_md.abfe_self_heal_cli --run-dir RUN/fep \\
        --decharge-lambdas ... --vdw-lambdas ...   # or rely on config defaults
    # modes:
    #   (default)   resubmit each unhealthy window via bsub
    #   --local     run pmemd.cuda in-place
    #   --dry-run   diagnose + edit mdin only, no execution

Exit code is 0 iff every window of every leg ends complete.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Import siblings via relative import when run as part of the amber_md package
# (python -m amber_md.abfe_self_heal_cli), but fall back to absolute/top-level
# imports so the module is also importable standalone (e.g. unit tests that add
# this directory to sys.path and do `import abfe_self_heal_cli`).
try:
    from .abfe_self_heal import heal_leg, prod_is_complete
    from .logger import get_logger
except ImportError:  # pragma: no cover - standalone/test use
    from abfe_self_heal import heal_leg, prod_is_complete
    try:
        from logger import get_logger
    except ImportError:
        get_logger = None
try:
    from .version import lsf_banner, version_string
except ImportError:  # pragma: no cover
    try:
        from version import lsf_banner, version_string
    except ImportError:
        def lsf_banner(ctx=""): return "# (amber_md version banner unavailable)\n"
        def version_string(): return "amber_md (version unknown)"
    def version_string(): return "amber_md (version unknown)"

try:
    log = get_logger(__name__)
except Exception:  # pragma: no cover - standalone use
    import logging
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger(__name__)


# The four ABFE legs. Solvent legs have no Boresch restraint but are otherwise
# identical in structure, so the same heal logic applies to all of them.
ABFE_LEGS = ("complex_decharge", "complex_vdw",
             "solvent_decharge", "solvent_vdw")


_PHYS_BLOWUP_CLASSES = {
    "blowup_temperature", "shake_failure", "box_drift",
    "softcore_com", "unknown", "missing_no_error",
}


def _clean_start_coord(wd: Path) -> str:
    """Cleanest available pre-production coordinate source."""
    return "min.rst" if (Path(wd) / "min.rst").exists() else "system.inpcrd"


def purge_poisoned_restarts(window_dir: Path, full_redo: bool) -> list:
    """v2.5.15: remove artifacts that would let pmemd resume (irest=1) from a
    DETONATED state. full_redo purges eq/dens restarts too (physics blow-up);
    otherwise only the incomplete prod.* (transient / external kill)."""
    wd = Path(window_dir)
    victims = ["prod.out", "prod.rst", "prod.nc"]
    if full_redo:
        victims += ["eq.rst", "eq.out", "eq.nc", "dens.rst", "dens.out", "dens.nc"]
    removed = []
    for name in victims:
        p = wd / name
        try:
            if p.exists():
                p.unlink(); removed.append(name)
        except OSError:
            pass
    return removed


def recovery_stages(window_dir: Path, failure_class: str = "unknown"):
    """v2.5.15: ordered (stage, mdin, start_coord, out, rst) for a recovery
    rerun, plus a full_redo flag.

    Physics blow-up  -> rebuild dens->eq->prod from a CLEAN coord source so the
                        remediation edits are exercised (the OLD code resumed
                        prod from the blown-up eq.rst and re-detonated).
    Transient/kill   -> resume prod from a known-good eq.rst if present.
    """
    wd = Path(window_dir)
    full_redo = failure_class in _PHYS_BLOWUP_CLASSES
    if full_redo:
        start = _clean_start_coord(wd)
        if (wd / "dens.in").exists():
            seq = [("dens", "dens.in", start, "dens.out", "dens.rst"),
                   ("eq", "eq.in", "dens.rst", "eq.out", "eq.rst"),
                   ("prod", "prod.in", "eq.rst", "prod.out", "prod.rst")]
        else:
            seq = [("eq", "eq.in", start, "eq.out", "eq.rst"),
                   ("prod", "prod.in", "eq.rst", "prod.out", "prod.rst")]
        return seq, full_redo
    if (wd / "eq.rst").exists():
        return [("prod", "prod.in", "eq.rst", "prod.out", "prod.rst")], full_redo
    start = _clean_start_coord(wd)
    return [("eq", "eq.in", start, "eq.out", "eq.rst"),
            ("prod", "prod.in", "eq.rst", "prod.out", "prod.rst")], full_redo


def _stages_in_to_run(window_dir: Path):
    """Back-compat shim: classify and return only the stage tuples."""
    try:
        from .abfe_self_heal import classify_failure
    except ImportError:  # pragma: no cover
        from abfe_self_heal import classify_failure
    fc = classify_failure(Path(window_dir)).failure_class
    seq, _ = recovery_stages(window_dir, fc)
    return seq


def make_local_rerun():
    """rerun callback that runs pmemd.cuda in-place for the failed stages."""
    try:
        from .abfe_self_heal import classify_failure
    except ImportError:  # pragma: no cover
        from abfe_self_heal import classify_failure

    def rerun(window_dir: Path, attempt: int) -> int:
        wd = Path(window_dir)
        fc = classify_failure(wd).failure_class
        seq, full_redo = recovery_stages(wd, fc)
        removed = purge_poisoned_restarts(wd, full_redo)
        log.info("  [local rerun] %s class=%s full_redo=%s purged=%s",
                 wd.name, fc, full_redo, removed)
        for stage, mdin, src, out, rst in seq:
            cmd = ["pmemd.cuda", "-O", "-i", mdin, "-p", "system.prmtop",
                   "-c", src, "-o", out, "-r", rst, "-x", f"{stage}.nc"]
            log.info("  [local rerun] %s attempt %d: %s", wd.name, attempt,
                     " ".join(cmd))
            rc = subprocess.run(cmd, cwd=wd).returncode
            if rc != 0:
                return rc
        return 0
    return rerun


_HOST_RE = re.compile(r"executed on host\(s\)\s+<([^>]+)>")


def _failed_hosts(window_dir: Path):
    """Return the set of hosts that previously RAN this window (and, since we
    only call this for an unhealthy window, killed it). Parsed from the LSF
    job-summary stdout files ('Job was executed on host(s) <gpu-node-02>').
    Scans the window dir and its parent leg dir (array .out lives at leg level).
    LSF host tokens can look like '1*gpu-node-02'; we strip the 'N*' prefix."""
    wd = Path(window_dir)
    cands = list(wd.glob("*.out"))
    if wd.parent.exists():
        cands += list(wd.parent.glob("fep_*.out")) + list(wd.parent.glob("heal_*.out"))
    hosts = set()
    for f in cands:
        if f.name in ("min.out", "dens.out", "eq.out", "prod.out"):
            continue
        try:
            txt = f.read_text(errors="replace")
        except OSError:
            continue
        for m in _HOST_RE.finditer(txt):
            for tok in m.group(1).split():
                hosts.add(tok.split("*")[-1].strip())
    return {h for h in hosts if h}


def _select_resource(avoid: set):
    """Build an LSF -R select expression that excludes the given hosts, e.g.
    -R "select[hname!='gpu-node-01' && hname!='gpu-node-02']". Returns '' if the
    avoid-set is empty."""
    hosts = sorted(h for h in avoid if h)
    if not hosts:
        return ""
    expr = " && ".join(f"hname!='{h}'" for h in hosts)
    return f'#BSUB -R "select[{expr}]"'


def make_bsub_rerun(project: str, queue: str, walltime: str,
                    mem_mb: int = 8192, avoid_hosts=None,
                    host_avoidance: bool = True):
    """rerun callback that resubmits a single window via bsub and BLOCKS until
    it finishes (-K), so the heal loop can re-check prod.out synchronously.

    Host-avoidance: each retry excludes (a) any host in the persistent
    `avoid_hosts` blocklist, PLUS (b) every host that previously ran this
    window (parsed from its LSF .out), so a window killed on gpu-node-02 is
    never rescheduled there. Disable with host_avoidance=False.
    """
    base_avoid = set(avoid_hosts or ())

    def rerun(window_dir: Path, attempt: int) -> int:
        wd = Path(window_dir)
        avoid = set(base_avoid)
        if host_avoidance:
            avoid |= _failed_hosts(wd)
        script = wd / "run_heal.lsf"
        lines = ["#!/bin/bash",
                 lsf_banner(f"self-heal rerun {wd.name} (attempt {attempt})").rstrip("\n"),
                 f"#BSUB -q {queue}", f"#BSUB -P {project}",
                 f"#BSUB -J heal_{wd.name}", f"#BSUB -W {walltime}",
                 "#BSUB -n 1",
                 '#BSUB -R "span[hosts=1]"',
                 f'#BSUB -R "rusage[mem={mem_mb}]"',
                 '#BSUB -gpu "num=1"',
                 f"#BSUB -o heal_{wd.name}.%J.out",
                 f"#BSUB -e heal_{wd.name}.%J.err"]
        sel = _select_resource(avoid)
        if sel:
            lines.append(sel)
            log.info("  [bsub rerun] %s avoiding host(s): %s",
                     wd.name, ", ".join(sorted(avoid)))
        lines += ["set -euo pipefail", "module purge"]
        for m in ("amber/22",):
            lines.append(f"module load {m}")
        lines.append(f'cd "{wd}"')
        try:
            from .abfe_self_heal import classify_failure
        except ImportError:  # pragma: no cover
            from abfe_self_heal import classify_failure
        _fc = classify_failure(wd).failure_class
        _seq, _full_redo = recovery_stages(wd, _fc)
        _removed = purge_poisoned_restarts(wd, _full_redo)
        log.info("  [bsub rerun] %s class=%s full_redo=%s purged=%s",
                 wd.name, _fc, _full_redo, _removed)
        for stage, mdin, src, out, rst in _seq:
            lines.append(f"pmemd.cuda -O -i {mdin} -p system.prmtop -c {src} "
                         f"-o {out} -r {rst} -x {stage}.nc")
        script.write_text("\n".join(lines) + "\n")
        script.chmod(0o755)
        log.info("  [bsub rerun] %s attempt %d (blocking)", wd.name, attempt)
        # 'bsub < script' needs the script on stdin; with shell=False we open
        # the file and pass it as the child's stdin (-K blocks until done).
        with open(script) as fh:
            rc = subprocess.run(["bsub", "-K"], stdin=fh, cwd=wd).returncode
        return rc
    return rerun

def make_dry_rerun():
    def rerun(window_dir: Path, attempt: int) -> int:
        log.info("  [dry-run] %s attempt %d: edits applied, NOT executing.",
                 Path(window_dir).name, attempt)
        return 0
    return rerun


def _leg_lambdas(leg_name, args, cfg):
    if leg_name.endswith("decharge"):
        return (args.decharge_lambdas
                or getattr(cfg, "decharge_lambdas", cfg.lambdas))
    return args.vdw_lambdas or getattr(cfg, "vdw_lambdas", cfg.lambdas)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Self-heal ABFE windows -> prod.out")
    ap.add_argument("--run-dir", required=True, type=Path,
                    help="The 'fep' directory containing the four leg dirs.")
    ap.add_argument("--legs", nargs="+", default=list(ABFE_LEGS),
                    help="Subset of legs to heal (default: all four).")
    ap.add_argument("--decharge-lambdas", type=float, nargs="+", default=None)
    ap.add_argument("--vdw-lambdas", type=float, nargs="+", default=None)
    ap.add_argument("--max-attempts", type=int, default=3)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--local", action="store_true",
                      help="Run pmemd.cuda in-place instead of bsub.")
    mode.add_argument("--dry-run", action="store_true",
                      help="Diagnose + edit mdin only; do not execute.")
    # bsub resources (used unless --local/--dry-run)
    ap.add_argument("--project", default="your-project")
    ap.add_argument("--queue", default="gpu")
    ap.add_argument("--walltime", default="24:00")
    ap.add_argument("--mem-mb", type=int, default=8192)
    # --- v2.5.8: true host-avoidance for genuine node failures ---
    ap.add_argument("--avoid-hosts", nargs="+", default=[],
                    metavar="HOST",
                    help="Persistent blocklist of GPU hosts to NEVER schedule "
                         "reruns on (e.g. --avoid-hosts gpu-node-01 gpu-node-02). "
                         "Added to the per-window set of hosts that already "
                         "killed the window.")
    ap.add_argument("--no-host-avoidance", action="store_true",
                    help="Disable automatic exclusion of hosts that previously "
                         "ran/killed a window (the --avoid-hosts blocklist is "
                         "still honored).")
    args = ap.parse_args(argv)

    try:
        from .config import FEPConfig
        cfg = FEPConfig()
    except Exception:
        class _C:  # minimal fallback
            lambdas = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                       0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
            decharge_lambdas = lambdas
            vdw_lambdas = lambdas
        cfg = _C()

    if args.dry_run:
        rerun = make_dry_rerun()
    elif args.local:
        rerun = make_local_rerun()
    else:
        rerun = make_bsub_rerun(args.project, args.queue, args.walltime,
                                args.mem_mb,
                                avoid_hosts=args.avoid_hosts,
                                host_avoidance=not args.no_host_avoidance)

    all_ok = True
    grand = {}
    for leg in args.legs:
        leg_dir = args.run_dir / leg
        if not leg_dir.exists():
            log.warning("leg dir missing, skipping: %s", leg_dir)
            continue
        lams = _leg_lambdas(leg, args, cfg)
        log.info("=== Healing leg %s (%d windows) ===", leg, len(lams))
        reps = heal_leg(leg_dir, lams, rerun, args.max_attempts, logger=log)
        grand[leg] = reps
        for lam, rep in reps.items():
            status = "OK" if rep.completed else "INCOMPLETE"
            log.info("  %s lambda_%s: %s (attempts=%d)",
                     leg, lam, status, rep.attempts)
            if not rep.completed:
                all_ok = False

    # concise summary
    print("\n==== SELF-HEAL SUMMARY ====")
    for leg, reps in grand.items():
        n_ok = sum(1 for r in reps.values() if r.completed)
        print(f"  {leg:18s}: {n_ok}/{len(reps)} windows complete")
        for lam, r in reps.items():
            if not r.completed:
                last = r.history[-1] if r.history else {}
                print(f"      lambda_{lam} STILL FAILING "
                      f"(class={last.get('failure_class','?')})")
    print("===========================")
    return 0 if all_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())