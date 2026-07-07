# Clean-run guide — amber_md v2.5.74

This is the deploy-and-launch procedure for a **clean ABFE run** on the
**8-GPU the GPU queue queue for ~1 week**. v2.5.74 is tuned for that budget.

## 0. What changed vs the last run (why a clean run is required)
The finished v2.5.72 `complex_decharge` leg was numerically STABLE (0 NaN, 0
drift) but NOT converged: TI −59 vs BAR −81 kcal/mol, MBAR unsolvable
(6504 foreign energies > 0 = overlap collapse), forward/reverse halves
~16 kcal/mol apart. Root causes are now fixed:

| Fix (v2.5.73 + v2.5.74) | What |
|---|---|
| Production length | 2 ns → **10 ns/window** |
| decharge schedule | 21 → **30 windows** (charging well 0.72–0.93 now ≤0.0375 λ spacing) |
| dense vdw schedule | → **50 windows** (soft-core danger zone 0.55–0.90 ≤0.0125); `use_dense_vdw=True` by default |
| convergence_analysis.py | None-guard (no more `len(None)` crash at λ=1.0) |
| analyze LSF | now actually RUNS `convergence_analysis.py` → writes `convergence.csv` / `.png` |
| HREMD builder | `cfg.lambdas` → `_active_lambdas` (correct window count) |

**Schedules changed, so any partially-completed leg is incompatible** (MBAR needs
one consistent λ set). Start a fresh run directory.

## 1. Budget — fits the week with headroom
All four legs use the same stage schedules (selected by stage, not by
complex/solvent): every *decharge* leg = 30 windows, every *vdw* leg = 50.

| Leg | Windows | ns/win | Leg ns | Wall @8 GPU |
|---|---|---|---|---|
| complex_decharge | 30 | 10 | 300 | ~28 h |
| complex_vdw | 50 | 10 | 500 | ~47 h |
| solvent_decharge | 30 | 10 | 300 | ~7 h (small box ~4× faster) |
| solvent_vdw | 50 | 10 | 500 | ~12 h |
| **TOTAL** | **160** | | **1600 ns** | **~3.9 days** |

≈ 750 GPU-hours = **56 % of an 8-GPU week**, leaving **~3 days headroom** to
extend any window the analyzer flags. (Estimate assumes ~32 ns/day/window
complex throughput from the prior run; solvent ~4× faster.)

## 2. Deploy (do NOT overwrite the old install)
```bash
cd ~/Tools
unzip amber_md_workflow_v2.5.74.zip          # -> amber_md_workflow_v2.5.74/
# point PYTHONPATH at the NEW path (the generated LSF bakes this in):
export PYTHONPATH=$HOME/Tools/affiniscreen
python -c "import amber_md; print(amber_md.__version__)"   # must print 2.5.74
```

## 3. Kill old jobs + start a fresh run dir
```bash
bkill 0                       # kill ALL your pending/running jobs (or bkill per-id)
# Use a NEW run dir; do not reuse the old abfe_* dir (schedule mismatch):
#   .../Run_dir/run_v250/abfe_<newdate>/...
```

## 4. Launch
Launch via your usual entrypoint (GUI or `run_amber.py`) against the new
PYTHONPATH. No flags needed — 10 ns, 30/50 windows and dense vdw are the
v2.5.74 defaults. (Legacy grids remain available: `--decharge-lambdas`,
`--vdw-lambdas`, or `use_dense_vdw=False` for the 40-window grid.)

## 5. Reading the convergence output (auto-generated per leg)
When a leg's analyze job runs it now writes, in the leg dir:
* `summary.json`  — headline ΔG + estimator
* `convergence.csv` — TI / BAR / BAR_err / MBAR / MBAR_err at 20/40/60/80/100 %
* `convergence.png` — the same as a plot
* `convergence.log` — full run log incl. fwd/rev halves and the VERDICT line

**Acceptance criteria for a leg:**
1. **MBAR solves** (MBAR column is a number, not `n/a`). `n/a` ⇒ overlap still
   too weak ⇒ that leg needs more/denser λ.
2. **TI ≈ BAR ≈ MBAR** — agreement within ~1–2 kcal/mol. A large TI–BAR gap is
   the classic "not converged / bad overlap" signature.
3. **|dG_fwd − dG_rev| ≲ 1 kcal/mol** (Q2 in the log). HIGH ⇒ undersampled.
4. **Plateau drift** (Q1 last-step) small ⇒ dG has stopped moving with time.

The log prints a one-line VERDICT (`converged` vs `NOT demonstrably
converged`). Trust it.

## 6. Extending a single non-converged window (use the headroom)
If only a few windows fail the fwd/rev test (usually around the λ≈0.8 well),
extend just those rather than re-running the leg:
```bash
cd <leg_dir>/lambda_0.800
# resume production from the last restart for N more ns (e.g. another 10 ns):
pmemd.cuda -O -i prod.in -p system.prmtop -c prod.rst -ref prod.rst \
           -o prod.ext.out -r prod.ext.rst -x prod.ext.nc
# then concatenate / re-point and re-run the analyze job for the leg.
```
After extending, re-run the leg's `analyze_<leg>.lsf` to regenerate
`convergence.csv` and re-check the criteria.

## 7. If λ≈0.8 is STILL not converged after 10 ns
That points to slow sampling (e.g. pocket-water reorganization), not λ spacing.
Options, in order: (a) extend those windows to 20 ns; (b) enable HREMD for that
leg **only if** you can place one GPU per replica on a single host (see
docs/HREMD_NOTES.md — Amber HREMD is synchronous and needs ~1 GPU/replica);
(c) run that ligand in OpenFE (async REX) instead.

## 8. Sanity checklist before you walk away
- [ ] `amber_md.__version__ == 2.5.74`
- [ ] new (empty) run dir, old jobs killed
- [ ] first few windows reach production (check a `prod.out` shows NSTEP climbing)
- [ ] density ~1.0 g/cc, T ~298 K in early prod output
- [ ] you know where `convergence.csv` will land (each leg dir)
