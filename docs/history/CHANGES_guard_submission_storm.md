# v2.5.0 build final47: stop the guard's duplicate-submission storm

## Symptom
`qu | wc -l` = 889 GPU jobs. One full ABFE pipeline = 76 array windows + 4
analyze + 1 cycle_close = 81 jobs. 889 / 81 ~= 11 -> the guard had stacked ~11
IDENTICAL full submissions in the queue.

## Root cause (guard design flaw, owned)
abfe_resume_guard.sh v1/v2 ran `fep_driver --resume --submit` every POLL_MIN.
`--resume --submit` ALWAYS creates a new LSF array for incomplete windows; it is
NOT pending-aware. On a busy gpu queue where nothing had finished yet, every
15-min sweep resubmitted all 76 windows again -> ~11 sweeps -> 889 jobs, all
duplicates of the same windows.

## Fix -- abfe_resume_guard.sh v3 is SUBMISSION-AWARE
* New live_jobs(): before any resubmit, count THIS run's PEND/RUN jobs, using
  the LSF array IDs the driver already records in <wd>/fep/job_ids.json
  (fallback: my fep_* jobs by name).
* THE GATE: if live_jobs > 0 -> SKIP the sweep entirely ("already queued,
  waiting"). Only resubmit when the queue has genuinely DRAINED (preemption /
  maintenance) AND windows remain incomplete. This makes duplicate stacking
  impossible.
* Stall accounting only advances when the queue is empty (real-crash detection
  is now meaningful, not triggered by a slow queue).
* Keeps final46 behaviour: passes --ligand-file + --ligand-resname on resume;
  FIRST_SWEEP_DELAY_MIN (30) before the first sweep.

## Operational
Before relaunch, clear the stacked jobs:
    pkill -9 -f abfe_resume_guard
    bkill 0
Then relaunch with the v3 guard (FIRST_LAUNCH=1 LIG_RESNAME=LIG ...).

## Driver fixes from final46 (carried) still apply:
  persisted-resname reload + fail-fast 0-atom mask guard.

## Verified
* py_compile clean (no .py changed vs final46; driver fixes intact).
* wrapper structurally validated; live-jobs gate + job_ids.json read present.
