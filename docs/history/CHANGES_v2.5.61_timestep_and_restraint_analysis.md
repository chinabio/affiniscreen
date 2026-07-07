# v2.5.61 — timestep stability fix + restraint-leg analysis closure

**Release date:** 2026-06-23  
**Base:** v2.5.60  
**Type:** functional fix (all FEP legs) + analysis findings + version-track resync

## 1. Summary
1. Production timestep dt 2 fs -> 1 fs (all legs), fixed via the GUI launch path.
2. Restraint-leg analysis closed: contribution is the analytic Boresch correction (-11.45 kcal/mol).
3. Version track resynced (VERSION + README banner were stale at 2.5.51a) to 2.5.61.
4. docs/release_v2.5.61/ bundles the crash-investigation handoff and the apply_dt001_* patch scripts.

## 2. Functional fix — production dt 0.002 -> 0.001 ps
### 2.1 Evidence (lig_12944901, multi-seed A/B)
| Leg / window | dt=0.002 | dt=0.001 |
|---|---|---|
| complex_restraint lambda 0.10-0.40 | crashed 3/3 seeds | stable |
| complex_decharge lambda_0.175 | died | (stable) |
| complex_vdw | killed before full test | use dt=0.001 |

Reproduced across seeds -> a true integrator blow-up, not box drift/clash. Cost ~2x wall-clock per ns.

### 2.2 GUI decoupling bug also fixed
`gui/pages/0_Setup_and_Launch.py` never passed `--dt` (dt came from `FEPConfig.dt_ps`) yet computed `--nstlim-prod = complex_ns*1e6/2` with a hard-coded /2 (2 fs). Changing dt alone would silently halve simulated time. Fixed both.

### 2.3 Changes
| File | Change |
|---|---|
| amber_md/config.py | FEPConfig.dt_ps 0.002 -> 0.001 |
| amber_md/config.py | MDConfig.prod_dt_ps 0.002 -> 0.001 |
| gui/pages/0_Setup_and_Launch.py | ABFE+RBFE _cns*1e6/2 -> /1 |
| gui/pages/0_Setup_and_Launch.py | ABFE+RBFE _ens*1e6/2 -> /1 |
| gui/fep_common.py | dt widget default 0.002 -> 0.001 |

eq_dt_ps/heat_dt_ps were already 0.001 (unchanged).

### 2.4 Verification gate (before committing GPU)
```
grep -E 'dt=|nstlim=' <run>/lig_*/fep/complex_vdw/lambda_0.500/prod.in
# expect: dt=0.001  and  nstlim=10000000  (complex_ns=10)
```

## 3. Restraint-leg analysis — findings and verdict
### 3.1 Verdict
dG_restraint = -11.45 kcal/mol (analytic Boresch standard-state correction). Re-derived from on-disk lambda=1 &rst+geometry to 0.01; self-test reproduces Deng-Roux T4 (11.62). Simulated windows are NOT a usable TI source.

### 3.2 FINDING 1 — FEPAnalyzer cannot analyze the Option-A leg
A real run returns 0/N windows, dG=None. alchemlyb needs icfe=1 (clambda/DV/DL); Option-A is icfe=0 -> "No free energy section found ... clambda was None". _fallback_ti/_parse_dvdl also find nothing. The v2.5.32/2.5.33 generator/analyzer schism, confirmed in practice. (decharge/vdw/solvent are icfe=1 -> analyzer works there.)

### 3.3 FINDING 2 — not a sound numerical-TI source
dG_attach = integral_0^1 <RESTRAINT>/lambda dlambda has a 1/lambda endpoint singularity. tools/restraint_ti_crosscheck.py integrates 3 ways: naive trapz ~88, sqrt-sub ~63, ln+tail ~87, vs analytic 11.45. ~40% spread = singularity artifact; dominant mass is in the unsampled lambda<0.1 region. Use the analytic term.

### 3.4 Maintainer guard
_leg() does total += rest['dG_kcal_mol'] + bcorr. An Option-A leg must report dG_kcal_mol=0 (or omit) so a singular numerical value can never corrupt the total.

## 4. Recommended fix for the design gap (pick ONE)
- Option X (simplest): analytic-only restraint contribution; analyzer recognizes Option-A (icfe=0/no clambda) and returns dG=0 + dG_boresch_correction=<analytic>.
- Option Y: soft/nonlinear restraint schedule + densify lambda<0.1 + emit icfe=1-style dHdl. More work; rarely beats analytic.
Add a Boresch consistency assert either way.

## 5. Version track
Resynced VERSION + README banner (were 2.5.51a) to 2.5.61, with __init__/run_amber/Home. Verified with tools/check_version_sync.py (clean).

## 6. New tools and docs
| Path | Purpose |
|---|---|
| tools/restraint_ti_crosscheck.py | Demonstrates the restraint-leg 1/lambda singularity vs analytic -11.45. |
| docs/release_v2.5.61/FEP_dt002_crash_investigation_HANDOFF.md | Full crash-investigation handoff (dt=0.002 instability, evidence, decisions). |
| docs/release_v2.5.61/apply_dt001_GUI_abfe_fix.sh | Standalone patch script for the GUI ABFE/Amber launch path (already applied in this release). |
| docs/release_v2.5.61/apply_dt001_clean_run_fix.sh | Earlier standalone clean-run patch script (run_amber.py path). |

## 7. Upgrade notes
1. Unpack v2.5.61 over your Tools/ location (back up the old tree first).
2. python -c "import amber_md; print(amber_md.__version__)" -> 2.5.61.
3. RESTART the Streamlit GUI.
4. Launch a fresh ABFE/Amber run; run the 2.4 gate BEFORE committing GPU.
5. Restraint contribution = analytic -11.45; do NOT FEPAnalyze the Option-A leg.
