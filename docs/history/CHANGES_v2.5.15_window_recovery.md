CHANGES v2.5.15 (final71) -- WINDOW RECOVERY (primary), MBAR rank guard, dense vdw
===============================================================================

PRIMARY FIX: windows that detonate now actually RECOVER
-------------------------------------------------------
Symptom: complex_vdw windows failed repeatedly and self-heal "could not
recover" them.

Root cause (3 bugs in the self-heal rerun path, abfe_self_heal_cli.py):
  1. _stages_in_to_run re-ran eq ONLY if eq.rst was MISSING. A detonated
     window had already WRITTEN an eq.rst (pmemd ntwr intermediate restart)
     before blowing up, so eq was SKIPPED and production restarted from the
     BLOWN-UP coordinates -> re-detonation every attempt. The remediation
     edits to eq.in were written but never executed.
  2. The remediation plan edits dens.in (attempt >=2) but dens was never in
     the rerun sequence, so box-collapse fixes were never applied.
  3. Stale incomplete prod.out / poisoned prod.rst / eq.rst were never
     removed; with irest=1,ntx=5 pmemd happily CONTINUED from the garbage.

Fix:
  * recovery_stages(window_dir, failure_class): for a PHYSICS blow-up
    (blowup_temperature, shake_failure, box_drift, softcore_com, unknown,
    missing_no_error) rebuild dens->eq->prod from a CLEAN coord source
    (min.rst / system.inpcrd). Transient/external kills keep the cheap
    resume-from-good-eq.rst path.
  * purge_poisoned_restarts(): delete the poisoned restart/output files before
    the rerun so pmemd cannot resume from a detonated state.
  * make_local_rerun / make_bsub_rerun now classify, purge, then run the
    recovery stages (and log class/full_redo/purged).
  * GUARANTEED last-resort (remediation_plan, attempt >= 4): a very gentle
    fixed protocol (dt=0.0005, gamma_ln=10, Berendsen barostat taup=5,
    softened soft-core scalpha=0.7/scbeta=18) from clean coords, so a window
    completes rather than looping. Raise --max-attempts to >=4 to enable it.

Validated (pure-Python pmemd simulation; no GPU): the OLD path re-detonates
(prod from poisoned eq.rst); the NEW path completes (rebuild from min.rst);
external_kill resumes prod from a good eq.rst with no wasteful redo. See
tests/test_self_heal_recovery.py.

SECONDARY (analysis) -- only relevant AFTER windows complete
------------------------------------------------------------
  * MBAR rank guard (fep.py): SVD rank/cond of u_nk; MBAR fit only when
    rank==n_states and cond<=1e8, else flagged ill-posed (the -1477 artifact).
    BAR/TI unaffected. See tests/test_mbar_rank_guard.py.
  * vdw_lambdas_dense (config.py): opt-in 40-window schedule (use_dense_vdw)
    for better overlap/conditioning so MBAR can become well-posed. Legacy
    28-window schedule remains the default; switching requires a full vdw
    re-run.

Operational notes
-----------------
  * Recovery is the fix you need to get every window to complete. Run the
    self-heal CLI with --max-attempts 4 (or more) to engage the guaranteed
    last-resort. Analysis fixes then apply to the completed leg.
