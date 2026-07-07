"""amber_md/abfe_self_heal.py -- per-window self-healing for ABFE/FEP legs.

Goal (user request): for EVERY lambda window of EVERY leg
(complex_decharge, complex_vdw, solvent_decharge, solvent_vdw), guarantee a
COMPLETE prod.out is produced. If a window finished without a complete
prod.out, inspect its diagnostic files (.err, eq.out, prod.out), classify the
failure, adjust the offending parameters, and re-run -- escalating the fix on
each attempt -- until prod.out is complete or attempts are exhausted.

This complements (does not replace) the existing in-LSF `run_prod_with_restart`
(which only handles the GPU "box changed too much" case) and `incomplete_indices`
(which only detects+resubmits). Here we additionally *diagnose and remediate*.

Design notes
------------
* Pure-Python and dependency-free (stdlib only) so it is unit-testable against
  the real eq.out fixtures and can run on the analyzer/login node.
* The window's mdin files (min.in/dens.in/eq.in/prod.in) are edited IN PLACE
  via small, well-scoped namelist substitutions. A backup `<name>.in.orig` is
  kept the first time a file is edited.
* `rerun_cmd` is pluggable so the same logic works for: a local pmemd.cuda
  re-run, an `bsub` resubmission of a single window, or a dry-run in tests.

The four legs are identical in structure, so one routine covers all of them.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

# --------------------------------------------------------------------------- #
# 1. Completion detection                                                      #
# --------------------------------------------------------------------------- #
# pmemd writes these only on a clean end-of-run. Presence of ANY => complete.
_PMEMD_DONE_MARKERS = (
    "Final Performance Info",
    "Master Total wall time",
    "wallclock() was called",
    "A V E R A G E S   O V E R",   # AVERAGES block is emitted at normal end
)


def prod_is_complete(prod_out: Path) -> bool:
    """True iff prod.out exists, is non-empty, and bears a pmemd end marker."""
    p = Path(prod_out)
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        txt = p.read_text(errors="replace")
    except OSError:
        return False
    return any(m in txt for m in _PMEMD_DONE_MARKERS)


# --------------------------------------------------------------------------- #
# 2. Failure classification                                                    #
# --------------------------------------------------------------------------- #
# Ordered most-specific -> most-generic. Each entry:
#   (class_name, [regexes], "which files to scan")
_SIGNATURES = [
    # soft-core / integrator blow-up (the complex_vdw lambda~0.5 case)
    ("blowup_temperature", [
        r"vlimit exceeded",
        r"TEMP\(K\)\s*=\s*(?:[1-9]\d{3,}|\d{5,})",   # T >= 1000 K
        r"Etot\s*=\s*\*{4,}",                         # energy overflow ****
        r"NaN",
    ], ("prod.out", "eq.out", "err")),
    # NPT box collapse / GPU grid invalidation
    ("box_drift", [
        r"box (?:dimensions|size) (?:have )?changed too much",
        r"changed too much from their initial",
        r"Periodic box dimensions have changed too much",
    ], ("prod.out", "err")),
    # fully-softcore COM-motion abort (solvent vdw high lambda)
    ("softcore_com", [
        r"completely softcore and skipped for C\.O\.M\.",
        r"Molecule \d+ is completely softcore",
    ], ("prod.out", "eq.out", "err")),
    # SHAKE failure -> usually too-large dt or bad start geometry
    ("shake_failure", [
        r"SHAKE",
        r"Coordinate resetting cannot be accomplished",
    ], ("prod.out", "eq.out", "err")),
    # external scheduler/node kill: pmemd ran fine then the JOB was terminated
    # (LSF "Exited with exit code 141", TERM_*, "Terminated at"). Distinct from
    # OOM: no memory signature, and pmemd wrote no error. Scan LSF .out summary,
    # the leg .err, and prod/eq mdout. (Driven by real gpu-node-02 exit-141.)
    ("external_kill", [
        r"Exited with exit code 1(?:41|43|37)",  # 128+SIG{TERM,KILL,INT}
        r"\bTERM_(?:RUNLIMIT|OWNER|PREEMPT|FORCE|OTHER)\b",
        r"Terminated at .* 20\d\d",               # LSF summary 'Terminated at'
    ], ("lsf",)),
    # host OOM / scheduler kill (truncated, no marker, no pmemd error)
    ("oom_or_killed", [
        r"Killed",
        r"out-of-memory|OOM|oom-kill",
        r"TERM_MEMLIMIT|TERM_RUNLIMIT|exceeded .* memory",
        r"DUE TO TIME LIMIT|CANCELLED",
    ], ("err",)),
    # generic input error
    ("input_error", [
        r"Input errors",
        r"unknown (?:flag|namelist)",
        r"Error opening",
    ], ("prod.out", "eq.out", "err")),
]


@dataclass
class Diagnosis:
    failure_class: str
    evidence: str = ""
    scanned: tuple = ()


def _read(p: Path) -> str:
    try:
        return Path(p).read_text(errors="replace")
    except OSError:
        return ""


def _err_files(window_dir: Path):
    """LSF writes fep_<leg>.<idx>.<job>.err at the LEG level, plus any *.err
    a window may have. Scan both the window dir and its parent (leg dir)."""
    wd = Path(window_dir)
    found = list(wd.glob("*.err"))
    if wd.parent.exists():
        found += list(wd.parent.glob("*.err"))
    return found
    return found


def _lsf_out_files(window_dir: Path):
    """LSF job-summary stdout files (fep_<leg>.<idx>.<job>.out, heal_*.out)
    carry the scheduler verdict -- 'Exited with exit code N', 'Terminated at',
    TERM_* -- which never appears in pmemd's .err. Scan window dir + leg dir."""
    wd = Path(window_dir)
    found = list(wd.glob("*.out"))
    if wd.parent.exists():
        found += list(wd.parent.glob("fep_*.out")) + list(wd.parent.glob("heal_*.out"))
    # exclude pmemd mdout files (min/dens/eq/prod.out) -- those aren't LSF logs
    return [f for f in found
            if f.name not in ("min.out", "dens.out", "eq.out", "prod.out")]

def classify_failure(window_dir: Path) -> Diagnosis:
    """Inspect a window's diagnostic files and return the best-matching class.

    Returns failure_class='unknown' if nothing matched (still actionable: the
    remediation policy applies a conservative generic stabilisation)."""
    wd = Path(window_dir)
    blobs = {
        "prod.out": _read(wd / "prod.out"),
        "eq.out": _read(wd / "eq.out"),
        "err": "\n".join(_read(f) for f in _err_files(wd)),
        "lsf": "\n".join(_read(f) for f in _lsf_out_files(wd)),
    }
    for cls, patterns, files in _SIGNATURES:
        for fkey in files:
            text = blobs.get(fkey, "")
            if not text:
                continue
            for pat in patterns:
                m = re.search(pat, text, re.IGNORECASE)
                if m:
                    line = m.group(0).strip()
                    return Diagnosis(cls, evidence=f"[{fkey}] {line[:160]}",
                                     scanned=tuple(files))
    # Nothing matched. If prod.out is simply missing/truncated with no error,
    # treat as a transient interruption (re-run as-is first).
    if not (wd / "prod.out").exists():
        return Diagnosis("missing_no_error", evidence="prod.out absent; no error signature")
    return Diagnosis("unknown", evidence="prod.out present but incomplete; no signature")


# --------------------------------------------------------------------------- #
# 3. Remediation policy                                                        #
# --------------------------------------------------------------------------- #
# Each policy returns a dict of {target_in_file: {namelist_var: new_value}}.
# `attempt` (1-based) escalates the strength of the fix.
def remediation_plan(failure_class: str, attempt: int) -> dict:
    """Map (failure_class, attempt) -> per-file namelist edits.

    Escalation philosophy: gentler integrator + stronger coupling + softer
    restraint, intensifying with each attempt. These mirror the final60 eq
    fix and Amber's standard recovery guidance.
    """
    a = max(1, attempt)
    # v2.5.15 GUARANTEED last-resort: on a high attempt index, apply the
    # gentlest fixed protocol regardless of class so a window completes
    # rather than looping. The CLI rebuilds from clean coords (recovery
    # path), so these edits are actually exercised.
    if attempt >= 4:
        return {
            "dens.in": {"dt": 0.0005, "nstlim": 1000000, "barostat": 1,
                         "taup": 5.0, "gamma_ln": 10.0},
            "eq.in":   {"dt": 0.0005, "gamma_ln": 10.0, "barostat": 1,
                         "taup": 5.0, "tempi": "${TEMP0}",
                         "scalpha": 0.7, "scbeta": 18.0},
            "prod.in": {"dt": 0.0005, "gamma_ln": 10.0,
                         "scalpha": 0.7, "scbeta": 18.0},
        }

    if failure_class in ("blowup_temperature", "shake_failure"):
        # shrink timestep, strengthen thermostat, gentler barostat on eq+prod
        dt = {1: 0.001, 2: 0.0005, 3: 0.0005}.get(a, 0.0005)
        gamma = {1: 5.0, 2: 10.0, 3: 10.0}.get(a, 10.0)
        edits = {
            "eq.in":   {"dt": dt, "gamma_ln": gamma, "tempi": "${TEMP0}",
                        "barostat": 1, "taup": 5.0},
            "prod.in": {"dt": dt, "gamma_ln": gamma},
        }
        if a >= 2:
            # also lengthen density pre-equilibration to settle the box first
            edits["dens.in"] = {"nstlim": 500000, "barostat": 1, "taup": 5.0}
        if a >= 3:
            # last resort: soften the soft-core core so the vdW derivative
            # singularity is smoothed (scbeta up, scalpha up)
            edits["eq.in"].update({"scalpha": 0.7, "scbeta": 18.0})
            edits["prod.in"].update({"scalpha": 0.7, "scbeta": 18.0})
        return edits

    if failure_class == "box_drift":
        # longer + gentler density equilibration; Berendsen barostat on prod
        return {
            "dens.in": {"nstlim": 500000 * a, "barostat": 1, "taup": 5.0},
            "prod.in": {"barostat": 1, "taup": 5.0},
        }

    if failure_class == "softcore_com":
        # disable COM-motion removal for fully-decoupled molecules
        return {"eq.in": {"nscm": 0}, "prod.in": {"nscm": 0}}

    if failure_class == "external_kill":
        # The job was terminated by the scheduler/node (e.g. LSF "Exited with
        # exit code 141", TERM_*, node/GPU event) -- NOT a physics problem.
        # Correct action is a clean resubmit with NO mdin changes. Optionally
        # nudge the scheduler to avoid the offending host on later attempts.
        plan = {"__resubmit__": {"reason": "external_kill", "edit_mdin": False}}
        if a >= 2:
            plan["__resources__"] = {"avoid_last_host": True,
                                     "walltime_bump": True}
        return plan

    if failure_class == "oom_or_killed":
        # not a parameter problem -> signal "resubmit with more resources"
        # (handled by caller; no mdin edit). Returned as a sentinel.
        return {"__resources__": {"fep_mem_mb": 8192 * (a + 1),
                                  "walltime_bump": True}}
    if failure_class == "input_error":
        return {}   # cannot auto-fix an input error safely; escalate to human

    if failure_class == "missing_no_error":
        # transient: first re-run as-is; only stabilise from attempt 2.
        if a == 1:
            return {}
        return remediation_plan("blowup_temperature", a - 1)

    # unknown: apply a conservative one-notch stabilisation
    return {"eq.in": {"dt": 0.001, "gamma_ln": 5.0},
            "prod.in": {"dt": 0.001, "gamma_ln": 5.0}}


# --------------------------------------------------------------------------- #
# 4. mdin editing (scoped namelist substitution)                              #
# --------------------------------------------------------------------------- #
def apply_mdin_edits(in_path: Path, var_values: dict) -> list:
    """Edit Amber namelist variables in place. Returns list of (var, old, new).

    Substitutes `var=...,` tokens within the &cntrl namelist. Adds the var
    before the terminating `/` if absent. Keeps a one-time `.orig` backup.
    Values that look like shell refs (e.g. '${TEMP0}') are written verbatim.
    """
    p = Path(in_path)
    if not p.exists():
        return []
    text = p.read_text()
    backup = p.with_suffix(p.suffix + ".orig")
    if not backup.exists():
        shutil.copy2(p, backup)

    changes = []
    for var, val in var_values.items():
        if var.startswith("__"):
            continue
        valstr = val if isinstance(val, str) else (
            f"{val:g}" if isinstance(val, float) else str(val))
        # match  var = <something> ,   (case-insensitive, optional spaces)
        pat = re.compile(rf"(\b{re.escape(var)}\s*=\s*)([^,\n/]+)(,?)",
                         re.IGNORECASE)
        m = pat.search(text)
        if m:
            old = m.group(2).strip()
            if old != valstr:
                text = pat.sub(lambda mm: f"{mm.group(1)}{valstr},", text, count=1)
                changes.append((var, old, valstr))
        else:
            # insert before the first namelist terminator '/'
            term = re.search(r"^\s*/\s*$", text, re.MULTILINE)
            insertion = f"  {var}={valstr},\n"
            if term:
                text = text[:term.start()] + insertion + text[term.start():]
            else:
                text = text.rstrip() + "\n" + insertion
            changes.append((var, "<absent>", valstr))
    if changes:
        p.write_text(text)
    return changes


# --------------------------------------------------------------------------- #
# 5. The self-heal driver                                                      #
# --------------------------------------------------------------------------- #
@dataclass
class HealReport:
    window_dir: str
    completed: bool
    attempts: int
    history: list = field(default_factory=list)   # list of dicts


def heal_window(window_dir: Path,
                rerun_cmd: Callable[[Path, int], int],
                max_attempts: int = 3,
                logger=None) -> HealReport:
    """Ensure `window_dir/prod.out` is complete, diagnosing+remediating+rerunning.

    Parameters
    ----------
    window_dir : the lambda_X.XXX directory.
    rerun_cmd  : callable(window_dir, attempt) -> int exit code. Caller decides
                 HOW to re-run (local pmemd, bsub single-window, or a test stub).
                 It must (re)generate prod.out.
    max_attempts : remediation attempts before giving up.
    """
    wd = Path(window_dir)

    def _log(msg):
        if logger:
            logger.info(msg)

    report = HealReport(window_dir=str(wd), completed=False, attempts=0)

    if prod_is_complete(wd / "prod.out"):
        report.completed = True
        _log(f"[self-heal] {wd.name}: already complete; nothing to do.")
        return report

    for attempt in range(1, max_attempts + 1):
        report.attempts = attempt
        diag = classify_failure(wd)
        plan = remediation_plan(diag.failure_class, attempt)

        applied = {}
        needs_resources = "__resources__" in plan
        for fname, edits in plan.items():
            if fname.startswith("__"):
                continue
            ch = apply_mdin_edits(wd / fname, edits)
            if ch:
                applied[fname] = ch

        step = {
            "attempt": attempt,
            "failure_class": diag.failure_class,
            "evidence": diag.evidence,
            "edits": applied,
            "resource_bump": plan.get("__resources__") if needs_resources else None,
        }
        _log(f"[self-heal] {wd.name}: attempt {attempt} class={diag.failure_class} "
             f"edits={ {k: len(v) for k, v in applied.items()} } "
             f"{'(+resource bump)' if needs_resources else ''}")

        if diag.failure_class == "input_error":
            step["action"] = "ABORT: input error needs human review"
            report.history.append(step)
            _log(f"[self-heal] {wd.name}: input error -- not auto-fixable.")
            break

        rc = rerun_cmd(wd, attempt)
        step["rerun_rc"] = rc
        complete = prod_is_complete(wd / "prod.out")
        step["complete_after"] = complete
        report.history.append(step)

        if complete:
            report.completed = True
            _log(f"[self-heal] {wd.name}: COMPLETE after attempt {attempt}.")
            break

    if not report.completed:
        _log(f"[self-heal] {wd.name}: still incomplete after {report.attempts} "
             f"attempt(s); leaving for human review.")
    return report


def heal_leg(leg_dir: Path,
             lambdas,
             rerun_cmd: Callable[[Path, int], int],
             max_attempts: int = 3,
             logger=None) -> dict:
    """Run heal_window over every lambda window of a leg. Returns
    {lambda_str: HealReport}. Applies uniformly to all four ABFE legs."""
    leg = Path(leg_dir)
    out = {}
    for lam in lambdas:
        wd = leg / f"lambda_{lam:.3f}"
        if not wd.exists():
            continue
        out[f"{lam:.3f}"] = heal_window(wd, rerun_cmd, max_attempts, logger)
    return out