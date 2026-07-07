# Convergence budget — the two levers (v2.5.75)

Two independent knobs control convergence. Diagnose which one a leg needs from
its `convergence.csv` / `convergence.log`:

| Symptom in convergence output | Failing lever | Fix |
|---|---|---|
| MBAR = `n/a` (won't solve); large TI-BAR gap | **lambda spacing** (overlap) | add windows in the steep region (`--decharge-lambdas` / `--vdw-lambdas`) |
| MBAR solves but |dG_fwd - dG_rev| HIGH; dG still drifting with time | **sampling time** | more ns/window (extend, see CLEAN_RUN) |

## Current defaults (v2.5.75)
* Production: **20 ns/window** (`nstlim_prod = 20_000_000` @ 1 fs)
* decharge: **30 windows** (charging well 0.72-0.93 <=0.0375 spacing)
* vdw: **50 windows** dense (soft-core danger zone 0.55-0.90 <=0.0125)
* Totals: 160 windows, 3200 ns, ~7.8 days wall @ 8 GPU, ~1500 GPU-h.

## Budget headroom
| Budget | GPU-h | This run uses | Headroom |
|---|---|---|---|
| 10 days | 1920 | 78% | ~2.2 d |
| 12 days | 2304 | 65% | ~4.2 d |
| 14 days | 2688 | 56% | ~6.2 d |

Use the headroom to extend ONLY the windows the analyzer flags (usually the
lambda~0.8 charging well), rather than re-running whole legs. See
docs/CLEAN_RUN_v2.5.74.md section 6 for the single-window extension recipe
(the procedure is identical for v2.5.75; only ns/window differs).

## If a leg still won't converge after 20 ns
The bottleneck is sampling QUALITY, not quantity. Options:
1. Extend the flagged windows to 40 ns (you have the headroom).
2. HREMD for that leg -- only if one GPU per replica fits on a single host
   (Amber HREMD is synchronous; see docs/HREMD_NOTES.md).
3. Run that ligand in OpenFE (async replica exchange).
