# v2.5.0 build final39: TI dV/dl endpoint-singularity guard

## Symptom
ABFE smoke test closed the cycle but dG_bind = +95.29 kcal/mol, UNTRUSTED.
estimator_spread ~1500-1800 kcal/mol on the complex legs; all 4 legs flagged
non-converged. Headline estimator = TI.

## Root cause (from fep/complex_vdw/lambda_1.000/prod.out)
At the fully-decoupled vdw endpoint (lambda=1.000) pmemd emitted an
endpoint singularity:
    Energy at 0.0000 = ****************      (overflow -> non-finite dV/dl)
    Etot = 0.0000 / 446.45 ...              (ghost end-state frames)
These blown-up dV/dl frames flowed straight into the TI estimator. The
existing _UNK_SANE_MAX_KT filter only cleaned the MBAR u_nk matrix, NOT the
TI dHdl stream -- and the headline was TI. One bad endpoint window therefore
poisoned the whole TI integral (-> +95 dG_bind).

## Fix (amber_md/fep.py)
* New FEPAnalyzer._sanitize_dHdl(dHdl): per-lambda removal of
    (a) non-finite dV/dl frames, and
    (b) robust outliers ( |x-median| > 25 * 1.4826 * MAD ), with an absolute
        cap |dV/dl| <= _DVDL_SANE_MAX_KCAL (5000 kcal/mol).
  Mirrors _sanitize_u_nk but for the gradient stream that TI actually uses.
* Applied right after dHdl = pd.concat(...), before statistical_inefficiency
  and TI(). Logs how many frames were dropped.
* New constants: _DVDL_SANE_MAX_KCAL=5e3, _DVDL_OUTLIER_MAD=25.

## Verified (synthetic dHdl reproducing the exact pathology)
  4 endpoint frames (inf, 1e6, -8000, 450 kcal at lambda=1.000) -> all dropped.
  Per-lambda dV/dl restored to 2.07 / 7.99 / 15.13 / 19.95.
  TI integral: BEFORE (poisoned) ~945  ->  AFTER (clean) ~8.9 kcal/mol.
  All .py compile.

## NOTE - this guards the ANALYSIS. The deeper sampling cure for the real run
## is still: more vdw windows (12-16, dense near lambda->1) + longer production
## (3-5 ns/window), and consider NOT sampling exactly lambda=1.000. The guard
## ensures a single bad endpoint frame can no longer produce a +95 headline.
