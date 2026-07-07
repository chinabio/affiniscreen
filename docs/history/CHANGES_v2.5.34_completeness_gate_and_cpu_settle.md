# v2.5.34 - reporting completeness gate + checkpoint resume + CPU density-settle

Date: 2026-06-20

## Fixes
1. REPORTING: run()/_collect_dvdl() gate each window on _window_complete() (pmemd
   end-of-run marker). Truncated prod.out => MISSING. New keys: short_windows,
   leg_failed, dG_reliable. Incomplete + ill-posed MBAR + diverged BAR => LEG FAILED.
2. RECOVERY: resume from prod.rst; box-drift ladder GPU regrid -> CPU density-settle
   (~25 ps) -> resume GPU. Engine auto-detect mpirun pmemd.MPI -> pmemd -> sander.
   Context-aware core cap: 32 on GPU nodes, 48 on normalQ. AMBERMD_CPU_NP overrides;
   else LSB_DJOB_NUMPROC. Loud [CPU_FALLBACK] logging.
   NOTE: the resume builder (build_lsf_array_resume) reuses build_lsf_array's body
   verbatim, so resuming the 6 missing windows automatically inherits BOTH the new
   harness and the new -n.
3. GPU FEP array job reserves CPU cores via NEW config field HPCConfig.fep_gpu_cores
   (default 16, was hardcoded 1). This builder requests its GPU via the explicit
   -gpu "num=" line, so -n is purely CPU cores -- no GPU-count inflation. GPU count
   unchanged.

## Clean version bump (2.5.34)
Read-then-write bump. check_version_sync.py passes across all 5 touchpoints.

## Deferred
Dense top-of-schedule deferred; rerun missing windows with the resume harness first.
