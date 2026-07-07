# v2.5.0 build final26: resume + auto-requeue wired into the GPU-queue path

## Goal
Make a requeued / crashed GPU job skip already-completed build + 5 ns MD and
resume from production.nc in the SAME workdir.

## What changed (_submit_to_lsf)
1. Inner command made EXPLICIT + deterministic. We strip any pre-existing
   --submit / --resume / --no-resume tokens and re-append exactly one canonical
   pair: `--submit local` plus `--resume` (or `--no-resume` if the user asked).
   Because the BSUB script does `cd <workdir>` and run_pipeline's resume guard
   skips build+MD when solvated.prmtop + production.nc exist, a re-run lands in
   the same dir and resumes (scores only).

2. LSF auto-requeue headers (default ON, --no-requeue to disable):
       #BSUB -Q "all ~0"   # requeue on any non-zero exit (0 = success, no requeue)
       #BSUB -r            # rerunnable: requeue on host failure / preemption

3. Bounded requeue to prevent an infinite resubmit loop on a DETERMINISTIC
   failure. The script tracks attempts in <workdir>/.requeue_attempts:
     * attempts 1..max_requeue: non-zero exit -> LSF requeues, resumes
     * attempt > max_requeue:   exit 0 (NOT in the -Q set) + write
                                REQUEUE_GIVING_UP.txt marker -> no further retry
     * on SUCCESS (exit 0):     counter + marker are cleared (set -e ensures the
                                cleanup line is only reached on success)
   New flag: --max-requeue N (default 3).

## CLI flags added
  --resume / --no-resume   (resume guard; default resume)
  --requeue / --no-requeue (LSF auto-requeue; default on)
  --max-requeue N          (attempt cap; default 3)

## Verified
* Generated .lsf: passes `bash -n`; amber module still dropped; cd <workdir>.
* Token handling: A) default gpu -> inner has `--submit local --resume`, no
  `--submit gpu`, `-Q`/`-r` present; B) --no-resume propagates; C) --no-requeue
  omits -Q/-r but still resumes; D) idempotent (no duplicate --resume / --submit
  local even if argv already had them).
* Bounded-requeue (faithful bash simulation): deterministic failure requeues
  3x then gives up at attempt 4 (exit 0 + marker); transient-then-success clears
  the counter. No infinite loop.
