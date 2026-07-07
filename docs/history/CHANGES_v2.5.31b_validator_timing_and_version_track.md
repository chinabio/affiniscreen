# v2.5.31b -- validator timing/noise fix + version-track discipline

VALIDATOR
  * Two-phase driver gate: content checks BEFORE the run script is written, then
    the -ref check AFTER run_<leg>.lsf exists -- both before any bsub. Fixes the
    2.5.31a gap where the missing-ref check (the 2.5.27 bug) was inert.
  * check_run_script runs ONCE per leg, targets run_<leg>.lsf, ignores analyze_*/
    cycle_close_* scripts (kills the 16-44x warning spam).

VERSION TRACK (now enforced on EVERY bump)
  * Updated together: amber_md/__init__.py (truth), VERSION, run_amber.py,
    README.md banner, GUI Home.py "What's new" (was stale since 2.5.26).
  * tools/check_version_sync.py  -- CI guard, fails on drift.
  * tools/bump_version.py        -- one-shot consistent bump + reminders.

OPERATIONAL NOTE (lig_12944901 crash)
  The crashed run's .lsf was stamped v2.5.31a, but its heat.in carried the
  pre-2.5.29 inline &wt TEMP_0 card -> 'Invalid TYPE flag'. That means the
  COMPUTE ENV imported a stale amber_md, not this package. Reinstall this package
  into that env (README -> "Install / version sync") so the imported code matches
  the stamp. The gate could not catch it because the stale env also lacked the gate.
