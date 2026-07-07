# v2.5.0 build final40: real fix for the +95 / spread~1500 ABFE artifact

## Correction to final39
final39 guarded the TI dV/dl stream. Re-analysis on the cluster proved TI was
NEVER the problem (TI = -45.97/-49.70/+11.10/-1.26, all sane). The ~1500
"spread" came from the MBAR/BAR (u_nk) path. final39 is harmless but did not
address the cause. final40 fixes the actual bugs.

## Root cause (verified in code + user logs)
1. alchemlyb's amber parser restarts each window's `time` index at 0. After
   pd.concat across lambdas the whole-frame timeseries is NOT single/
   contiguous/sorted, so:
       statistical_inefficiency failed: 'Duplicate time values found ...'
   Decorrelation was skipped; MBAR/BAR then received malformed u_nk ->
       MBAR failed: SVD did not converge ; DLASCL illegal value
       BAR = -1551 / -1846  (garbage)
2. estimator_spread_kcal = |BAR - TI| used that GARBAGE BAR, so every leg was
   flagged estimator_converged=False -> cycle UNTRUSTED, dG_bind +95.

## Fixes
### fep.py
* FIX 1 (decorrelation): new _decorrelate_per_lambda() runs
  statistical_inefficiency PER lambda group (de-dup + sort each window) then
  concats. Eliminates the 'Duplicate time values' failure so MBAR/BAR get
  well-formed u_nk. (alchemlyb-recommended subsampling granularity.)
* FIX 2 (honest spread):
    - the estimator fit loop records exceptions in _failed_estimators.
    - an estimator whose |dG - TI| > _DIVERGENCE_KCAL (100) is treated as a
      NUMERICAL FAILURE (collapsed u_nk / SVD), marked diverged=True, and
      EXCLUDED from estimator_spread_kcal / estimator_converged.
    - only_ti_survived=True when BAR and MBAR both fail/diverge.
    - summary.json now carries failed_estimators + only_ti_survived.
* cycle-closer trust gate: a leg with only_ti_survived is added to
  mbar_failed_legs and EXCLUDED from trusted=true, with an explicit untrusted
  reason ("MBAR/BAR failed ... TI-only ...") instead of a phantom spread.
  ABFE_RESULT.json gains mbar_failed_legs.

### abfe_qc.py
* Distinguishes a TI-only MBAR FAILURE (now: "MBAR/BAR unavailable -- TI-only;
  typical of truncated/too-few-window runs") from a genuine TI/MBAR
  disagreement ("estimator spread X; add windows"). No more misleading
  "spread 1500 -> lengthen production".

## Verified
* Fix 1: raw concat is NOT single/contiguous/sorted (the exact raise
  condition); per-lambda split makes every window valid. (pandas check)
* Fix 2: with the user's numbers
    complex_decharge TI=-45.97, BAR=-1551 -> BAR excluded, only_ti=True,
      flagged MBAR_FAILED (NOT spread=1505).
    healthy leg (TI=-8.9,BAR=-9.1,MBAR=-8.8) -> spread 0.2, converged.
    real disagreement (TI=-8.9,BAR=-13) -> spread 4.1, NOT converged.
* All .py compile.

## HONEST EXPECTATION on re-analyzing the existing smoke test
The smoke test is truncated (e.g. lambda_0.750 = 38/100 frames) and has only
5 windows / 100 ps. After final40:
  - the false +95 / spread-1500 artifact disappears,
  - legs report TI-only with an MBAR-failed reason (honest),
  - the cycle remains UNTRUSTED -- for the RIGHT reason (MBAR unavailable on a
    truncated run), not a phantom number.
A trustworthy MBAR dG still requires the production run (more windows, longer
sampling, finished windows).
