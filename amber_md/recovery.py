"""v2.5.15: completeness gate + recovery driver + structured logging.

Used by fep_driver BEFORE any leg is analyzed. Guarantees we attempt to make
every window complete first, and records exactly what was done to recover each
window (per attempt: failure class, evidence, mdin edits, purged files, rerun
rc, completeness) into:

    fep/<leg>/recovery_log.json     (machine-readable, per window/attempt)
    fep/RECOVERY_REPORT.txt         (human-readable roll-up across legs)

so a later troubleshooting session can see the full remediation history.
"""
from __future__ import annotations
from pathlib import Path
import json, time
from dataclasses import asdict, is_dataclass

from .abfe_self_heal import heal_leg, prod_is_complete


def leg_completeness(leg_dir, lambdas):
    """Return (n_complete, n_total, missing) for a leg without re-running."""
    leg = Path(leg_dir)
    total, complete, missing = 0, 0, []
    for lam in lambdas:
        wd = leg / ("lambda_%.3f" % lam)
        if not wd.exists():
            continue
        total += 1
        if prod_is_complete(wd / "prod.out"):
            complete += 1
        else:
            missing.append("%.3f" % lam)
    return complete, total, missing


def _report_from_heal(reps):
    """Normalise {lam: HealReport} into a JSON-able dict."""
    out = {}
    for lam, r in reps.items():
        d = asdict(r) if is_dataclass(r) else dict(getattr(r, "__dict__", {}))
        out[lam] = d
    return out


def ensure_leg_complete(leg_dir, lambdas, rerun, max_attempts, logger=None):
    """Gate one leg: if any window is incomplete, run heal_leg and LOG every
    recovery action. Returns a dict with completeness + the full heal history.

    Writes fep/<leg>/recovery_log.json (always, when recovery is attempted).
    """
    leg = Path(leg_dir)
    pre_c, pre_n, pre_missing = leg_completeness(leg, lambdas)
    rec = {
        "leg": leg.name,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "pre": {"complete": pre_c, "total": pre_n, "missing": pre_missing},
        "recovery_attempted": False,
        "windows": {},
    }
    if pre_missing:
        if logger:
            logger.warning("  [recovery] %s: %d/%d windows complete; "
                           "attempting recovery of %s",
                           leg.name, pre_c, pre_n, pre_missing)
        rec["recovery_attempted"] = True
        reps = heal_leg(leg, lambdas, rerun, max_attempts, logger=logger)
        rec["windows"] = _report_from_heal(reps)
        post_c, post_n, post_missing = leg_completeness(leg, lambdas)
        rec["post"] = {"complete": post_c, "total": post_n,
                       "missing": post_missing}
        # log a concise per-window remediation trail
        if logger:
            for lam, wd in rec["windows"].items():
                for step in wd.get("history", []):
                    logger.info("  [recovery] %s lambda_%s attempt %s "
                                "class=%s edits=%s rc=%s complete=%s",
                                leg.name, lam, step.get("attempt"),
                                step.get("failure_class"),
                                step.get("edits"), step.get("rerun_rc"),
                                step.get("complete_after"))
            if post_missing:
                logger.error("  [recovery] %s: STILL INCOMPLETE after recovery: "
                             "%s -- leg will NOT be analyzed as trusted.",
                             leg.name, post_missing)
            else:
                logger.info("  [recovery] %s: ALL %d windows complete after "
                            "recovery.", leg.name, post_n)
    else:
        rec["post"] = {"complete": pre_c, "total": pre_n, "missing": []}
        if logger:
            logger.info("  [recovery] %s: all %d windows already complete; "
                        "no recovery needed.", leg.name, pre_n)
    try:
        (leg / "recovery_log.json").write_text(json.dumps(rec, indent=2,
                                                          default=str))
    except OSError:
        pass
    rec["complete"] = not rec["post"]["missing"]
    return rec


def write_recovery_report(fep_root, leg_records, logger=None):
    """Human-readable roll-up across all legs -> fep/RECOVERY_REPORT.txt."""
    lines = ["AMBER FEP recovery report",
             "generated %s" % time.strftime("%Y-%m-%dT%H:%M:%S"),
             "=" * 60]
    all_ok = True
    for rec in leg_records:
        post = rec.get("post", {})
        ok = not post.get("missing")
        all_ok = all_ok and ok
        lines.append("")
        lines.append("LEG: %s  ->  %s" % (rec["leg"],
                     "COMPLETE" if ok else "INCOMPLETE %s" % post.get("missing")))
        lines.append("  before: %s/%s complete" % (rec["pre"]["complete"],
                                                   rec["pre"]["total"]))
        lines.append("  after : %s/%s complete" % (post.get("complete"),
                                                   post.get("total")))
        if rec.get("recovery_attempted"):
            for lam, wd in rec.get("windows", {}).items():
                status = "OK" if wd.get("completed") else "FAILED"
                lines.append("    lambda_%s: %s (%d attempt(s))"
                             % (lam, status, wd.get("attempts", 0)))
                for step in wd.get("history", []):
                    lines.append("      attempt %s: class=%s edits=%s rc=%s "
                                 "complete=%s"
                                 % (step.get("attempt"),
                                    step.get("failure_class"),
                                    step.get("edits"),
                                    step.get("rerun_rc"),
                                    step.get("complete_after")))
                    if step.get("evidence"):
                        lines.append("        evidence: %s" % step["evidence"])
    lines.append("")
    lines.append("=" * 60)
    lines.append("OVERALL: %s" % ("ALL LEGS COMPLETE" if all_ok
                                  else "SOME WINDOWS STILL INCOMPLETE"))
    txt = "\n".join(lines) + "\n"
    try:
        (Path(fep_root) / "RECOVERY_REPORT.txt").write_text(txt)
    except OSError:
        pass
    if logger:
        logger.info("  [recovery] wrote %s/RECOVERY_REPORT.txt (overall: %s)",
                    fep_root, "COMPLETE" if all_ok else "INCOMPLETE")
    return all_ok
