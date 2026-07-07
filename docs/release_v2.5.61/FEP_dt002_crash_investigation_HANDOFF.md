# FEP Production Crash — Investigation & Rescue Handoff

**Date:** 2026-06-22 (updated 2026-06-23)
**System:** `run_v250/abfe_20260621_204011/lig_12944901` (ABFE, Amber 22 pmemd.cuda, L40S GPUs)
**Affected leg:** `fep/complex_restraint`
**Affected windows:** lambda_0.100, 0.144, 0.196, 0.256, 0.324 AND 0.400 (mid-lambda band)
**Status:** ROOT CAUSE FOUND & VERIFIED. 0.100-0.324 rescued+validated (5/5 PASS).
0.400 rescue submitted (rescue_restraint_0400.sh). Generator fix written. Boresch
physics verified trustworthy. See sections 8-11 (added 2026-06-23) for the latest.
---

## 1. TL;DR

- **Symptom:** mid-lambda `complex_restraint` production windows die with a one-step
  energy explosion (BOND -> ~1.2e6, then NaN), followed by pmemd's
  "Periodic box dimensions have changed too much" halt. The recovery ladder
  (regrid / CPU-settle) loops and never produces a usable trajectory.
- **Root cause (PROVEN):** the production timestep **`dt=0.002` is integrator-
  unstable for this system.** It sits on the stability edge, so a random Langevin
  kick tips one degree of freedom over the cliff in a single step — stochastic
  timing (observed 2 ps to 42 ps), seed-dependent.
- **Fix (VERIFIED):** run production at **`dt=0.001`**. Ran the full length at a
  steady 298 K under NVT, NPT-Berendsen, and NPT-MC barostats, 3 random seeds
  each (9 stable runs). dt=0.002 crashed 3/3.
- **Rescue in flight:** `rescue_dt001_array.sh` re-runs the 5 stuck windows at
  dt=0.001 (nstlim doubled to keep 5 ns), one GPU per window, from the existing
  healthy `eq.rst`. No re-equilibration needed.
- **Still open:** permanent generator fix in `amber_md.fep.FEPSetup` (set prod
  dt=0.001; audit other legs). Minor: Boresch force constants look ~0.14x
  canonical (accuracy, not stability).

---

## 2. Confirmed facts (the evidence base)

All established from on-disk data, not inference:

| Fact | How confirmed |
|---|---|
| eq.rst is HEALTHY (box 124.4 A, 1.482e6 A^3, angles 109.47) | `decode_eqrst_box.py` via scipy NetCDF read |
| eq ran 750 ps stable at density 1.012 | eq.out |
| prod starts hot (~299 K at NSTEP 1000) then crashes | prod.out |
| Crash is a 1-step BOND-term explosion (2600 -> 1.27e6) | prod.out NSTEP 1000->2000 |
| Crash timing is STOCHASTIC (2 ps one run, 42 ps another), same eq.rst/prod.in, only ig=-1 seed differs | two prod.out attempts |
| 5 failing windows ran on 5 DIFFERENT GPU nodes (gpu01,03,04,05,06) | prod.out Hostname lines |
| eq uses dt=0.001 + ntr=1 positional restraints + gamma_ln=5.0 | eq.in |
| prod uses dt=0.002, NO ntr, gamma_ln=2.0 | prod.in |
| No softcore/dummy atoms (no icfe/ifsc/timask/scmask) | grep of *.in |
| dt=0.002 NPT: 3/3 crash (halt, NaN, real density 1.1-1.6) | test arm A |
| dt=0.001 NPT/NVT/MC: 9/9 stable, full length, steady 298 K | test arms B,C,D,E |

### The "phantom record" gotcha (IMPORTANT for anyone parsing mdout)
Every Amber mdout here ends with a **spurious trailing record**: a third copy of
the final `NSTEP=...` line showing **TEMP ~= 0.75 K** and **Density ~= 0.0004-0.0009**.
This is NOT a dynamics frame — it is a benign output artifact present in ALL runs
(including healthy ones and the original eq.out).

**Do not parse mdout with `tail -1`.** Earlier test harnesses did, mistook this
phantom for the final state, and produced FALSE "box evaporated to vacuum"
verdicts for runs that were actually perfectly stable at 298 K. Always filter to
records with TEMP >= 100 K (the rescue script's awk parser does this).

---

## 3. Hypotheses tested and ELIMINATED (so we don't re-chase them)

| # | Hypothesis | Why it was ruled out |
|---|---|---|
| 1 | Box-drift / grid-tolerance false alarm | box genuinely changes; halt is real |
| 2 | Recovery ladder rewinding to origin restart | restart source irrelevant to crash |
| 3 | Barostat algorithm (Berendsen vs MC) | both stable at dt=0.001; both "fail" at dt=0.002 |
| 4 | Corrupt / vacuum eq.rst | decoded healthy: 1.482e6 A^3, 124.4 A edges |
| 5 | Angle restraint singularity at r4=180 | actual angles 70/96 deg, far from 180 |
| 6 | Hardware / GPU ECC fault | 5 failures on 5 different nodes |
| 7 | Boresch restraint reference mismatch | healthy windows are mismatched MORE than failing ones |
| 8 | Decoupled/dummy atoms without softcore | none exist (no icfe/ifsc) |
| 9 | Thermostat cooling system to 0 K | FALSE — was the phantom-record parsing artifact; system runs at 298 K |
| -> | **dt=0.002 integrator instability** | **CONFIRMED: dt=0.002 3/3 crash, dt=0.001 9/9 stable** |

Why dt=0.002 is the answer and everything else fits:
- eq always survived -> eq uses dt=0.001
- regrid "rescues" briefly worked -> they used dt 0.0005-0.001
- crash is 1-step BOND explosion -> classic timestep instability signature
- stochastic timing -> dt=0.002 on the stability edge; random kick tips it over
- node-independent, healthy box, both barostats -> consistent with integrator, not
  pressure/hardware/restraint/topology

---

## 4. The rescue (CURRENT APPROACH — Option A, pure GPU, no CPU fallback)

**Script:** `rescue_dt001_array.sh` (LSF job array, 1 window per element, 1 GPU each)
**Submitted as:** job 4701985 `rescue_dt001[1-5]`

Per-window actions:
1. Back up original `prod.in` -> `prod.in.dt002.bak` (once, non-destructive).
2. Write new `prod.in`: only change is `dt=0.002 -> 0.001`, `nstlim` doubled
   2.5e6 -> 5.0e6 (keeps total sampling at 5 ns). Barostat stays Berendsen
   (barostat=1), gamma_ln=2.0, skinnb=3.0 — unchanged from original prod.
3. `pmemd.cuda -O -i prod.in -p system.prmtop -c eq.rst -o prod.out -r prod.rst -x prod.nc`
4. Parse REAL last frame (TEMP>=100 K filter) and print PASS/FAIL.

Array index -> window map: 1=0.100, 2=0.144, 3=0.196, 4=0.256, 5=0.324

**Note on cores/CPU fallback:** this rescue is pure GPU (pmemd.cuda), so it uses
only 1 GPU + a few host cores (`-n 4`). It deliberately does NOT use the 16-core
CPU-settle fallback, because dt=0.001 is verified stable and there is nothing to
fall back from. If a window unexpectedly fails, inspect it manually rather than
relying on a fallback. (The permanent generator fix SHOULD retain the 16-core
CPU-settle fallback for general robustness across other systems — see Open Items.)

### Monitoring / collecting results
```bash
bjobs -A                                  # array summary
bjobs                                      # per-element + node
grep -H "VERDICT" rescue_dt001.*.out       # collect all 5 verdicts when done
```
Expected: all 5 -> `VERDICT: PASS`, NSTEP=5000000, TEMP ~298, Density ~1.01.

If a window halts before NSTEP=5000000: note its NSTEP and inspect that window;
do not assume the phantom record (TEMP~0.75) means failure.

---

## 5. Open items (permanent fixes — NOT yet done)

1. **Generator fix in `amber_md.fep.FEPSetup` (PRIMARY):**
   - Set production `dt=0.001` for the complex_restraint leg.
   - Keep total sim time constant by doubling nstlim.
   - **Audit the other legs** for the same dt=0.002 risk: `complex_decharge`,
     `complex_vdw`, `solvent_decharge`, `solvent_vdw` (legs seen under `fep/`).
   - RETAIN the 16-core CPU-settle fallback in the generated workflow for general
     robustness (this rescue skips it only because we already proved dt=0.001
     stable for THIS system).

2. **Boresch force constants (ACCURACY, low priority):**
   - boresch.RST uses rk2=1.440 (distance) and 14.400 (angles/dihedrals),
     ~0.14x the canonical Boresch values (~10 kcal/mol/A^2 distance,
     ~100 kcal/mol/rad^2 angular). Confirm against intended FEPSetup values;
     too-weak restraints can bias dG (will not crash anything).

3. **mdout phantom trailing record:**
   - Document/teach the parsing pitfall (TEMP~0.75 K / density~0.0004 trailing
     line). Any QC tooling must filter TEMP>=100 K before reading "final" state.

---

## 6. Key files (in cc-file-api / on disk)

| File | Purpose |
|---|---|
| `rescue_dt001_array.sh` | CURRENT rescue (parallel, 1 GPU/window) — the one running |
| `rescue_dt001_all_windows.sh` | serial version (superseded; kept for reference) |
| `decode_eqrst_box.py` | reads box/volume from NetCDF restart (proved eq.rst healthy) |
| `boresch_match_check.py` | compares restraint refs vs eq geometry (ruled out hyp #7) |
| `test_nvt_vs_npt_lambda0144.sh` | the decisive NVT/NPT multi-seed test |
| `test_dt_fix_lambda0144.sh` | the dt A/B/C multi-seed test (proved dt=0.001 fix) |

On-disk working dir per window:
`run_v250/abfe_20260621_204011/lig_12944901/fep/complex_restraint/lambda_<L>/`
Originals preserved as `prod.in.dt002.bak`; `orig_eq.rst` is the pre-ladder eq restart.

---

## 7. If you must restart the investigation from scratch

You don't need to. The answer is dt=0.001. To re-verify in one command set:
```bash
WD=.../lambda_0.144
# confirm dt=0.002 crashes and dt=0.001 survives (see test_dt_fix_lambda0144.sh)
# confirm eq.rst healthy:
python decode_eqrst_box.py "$WD/eq.rst"     # expect VOLUME ~1.48e6
# confirm the protocol gap:
grep -hE "dt=|gamma_ln|ntr=" "$WD/eq.in" "$WD/prod.in"
#   eq.in: dt=0.001, gamma_ln=5.0, ntr=1   |   prod.in: dt=0.002, gamma_ln=2.0, (no ntr)
```

---

# ===== UPDATES 2026-06-23 =====

## 8. Rescue outcome + the 6th window the audit caught

- **Rescue of the original 5 windows: ALL PASS.** lambda_0.100/0.144/0.196/0.256/
  0.324 ran the full 5 ns at dt=0.001. Spot-checked lambda_0.144: real final
  frame NSTEP=5,000,000, TEMP 297-298 K, Density 1.0174-1.0177, prod.nc = 889 MB
  of real frames. The phantom record (TEMP=0.81 K) is present and correctly
  ignored by the parser.
- **The audit found a 6th failed window: lambda_0.400** (prod.nc only 17 MB, still
  dt=0.002). It died the same way; the original "5 stuck windows" count was wrong.
  Rescue submitted via `rescue_restraint_0400.sh` (dt=0.001, nstlim=5e6).
- **Restraint leg window map (dt + status):**
  | windows | dt | status |
  |---|---|---|
  | 0.000,0.004,0.016,0.036,0.064 | 0.002 | complete (848 MB) |
  | 0.100,0.144,0.196,0.256,0.324 | 0.001 | RESCUED, validated |
  | 0.400 | 0.002->0.001 | rescuing now |
  | 0.500..1.000 | 0.002 | complete (848 MB) |
  Note: dt=0.002 survived at 0.000-0.064 and 0.500-1.000 but died at 0.100-0.400
  -> the strained mid-attachment band is where the 2 fs step tips over. Fully
  consistent with "dt=0.002 on the stability edge."

## 9. Blast-radius audit (tools/audit_dt_all_legs.sh)

ALL legs write production at dt=0.002 (shared cfg.dt_ps). Per-leg disk status:
- complex_decharge: 15/15 ran; lambda_0.175 only 3 MB (died) -- NOT rescued (user
  is doing a clean run; out of scope here).
- complex_restraint: 35/36 complete; 0.400 rescuing (see sec 8).
- complex_vdw: only 0.000, 0.050, 0.100 have any data; rest empty/none.
- solvent_decharge, solvent_vdw: none run.
The decharge/vdw/solvent legs are dt=0.002 and SHARE THE RISK, but per user
decision we fix the RESTRAINT leg first and revisit the rest after a clean run.

## 10. complex_vdw: KILLED, not crashed (do NOT assume dt is the cause)

IMPORTANT correction to an earlier assumption. The vdw windows that stopped early
(0.000 @562k, 0.050 @409k, 0.100 @4k steps) show:
  * healthy last frame (BOND ~2600, TEMP ~298, small smooth dV/dl = 0/2.3/9.0)
  * NO crash fingerprints: no "changed too much", no NaN, no BOND explosion,
    no timing footer.
  * user confirms they bkill'd the vdw jobs.
=> These were KILLED mid-run, NOT crashed. There is currently NO evidence the
   vdw leg suffers the dt instability. The variable stop-steps are just where the
   kill landed, not stochastic crashes.
vdw IS a true softcore alchemical leg (icfe=1, ifsc=1, scalpha=0.2, scbeta=50,
gti_*, mbar_states=44) -- a DIFFERENT machine from the restraint leg. The package
itself warns of a softcore endpoint problem near high lambda (config.py:236
"0.15..0.70 hit prod box drift -- the classic soft-core"; analysis header
"dV/dl -> -113 at high lambda"). User recalls vdw "always blows up at 0.600".
PLAN: do NOT guess. After the restraint leg is fixed and a clean run is done,
diagnose vdw from REAL crash data (a genuine 0.600 failure with NaN/halt/BOND
explosion), not from killed-job artifacts. Discriminator:
  * BOND/velocity explosion from healthy frame -> timestep cliff -> dt=0.001 helps
  * dV/dl diverging to -100s near 0.6 -> softcore singularity -> needs
    scalpha/scbeta or use_dense_vdw/schedule fix, NOT dt.
The vdw dt A/B harness (test_dt_fix_vdw_0100.sh) was built but SET ASIDE -- there
is no confirmed crash to reproduce. Do not run it against killed-job windows.

## 11. Boresch physics: VERIFIED TRUSTWORTHY (earlier "0.14x" alarm was FALSE)

The earlier worry that force constants were ~0.14x canonical (a possible bug) is
RETRACTED. The constants are correct -- they are lambda-scaled on purpose.

Code path (amber_md/fep.py _write_boresch_RST, ~584-639; amber_md/boresch.py):
  * base canonical constants: kr=10.0, kth=100.0, kph=100.0 (boresch.py:94-96)
  * writer scales by lambda: _scale=float(lam); kr*=lam; kth*=lam; kph*=lam
    (fep.py:614-617), called with lam=window-lambda (fep.py:746).
  * restraint leg therefore RAMPS the restraint ON from 0 (lambda=0) to full
    canonical (lambda=1). This is intended (Option A: lambda-scaled &rst, icfe=0).

ON-DISK VERIFICATION (boresch.RST rk2/rk3 values):
  | window | distance rk | angle/dih rk | = lambda x (10,100)? |
  |---|---|---|---|
  | 1.000 | 10.000 | 100.000 | yes |
  | 0.500 |  5.000 |  50.000 | yes |
  | 0.144 |  1.440 |  14.400 | yes |
  Clean linear scaling. Geometry sane: r0=5.135 A, angles 88.4/100.8 deg
  (away from 0/180 singularities).

ANALYTIC CORRECTION (lambda_1.000/boresch_correction.txt): -11.447503 kcal/mol.
  * lambda=1 endpoint carries FULL canonical k (10/100/100), which is exactly
    what the analytic Deng-Roux formula uses -> simulated endpoint and analytic
    term are CONSISTENT.
  * -11.45 sits in the expected 11-12 kcal/mol band (package self-test:
    dG_release 11.63 vs Deng-Roux T4 paper 11.62). Sign correct (ADD = -release).
=> The restraint dG (TI over dV/dl on the scaled &rst + the -11.45 analytic
   standard-state correction) is SELF-CONSISTENT and trustworthy.

### 11a. LATENT BUG / TRAP (file for maintainer; does NOT affect current run)
The writer code (v2.5.32 lambda-scaling, STILL EXECUTING) contradicts its own
v2.5.33 comment (fep.py:611-613) which claims "FIXED Boresch force constants for
ALL stages... lambda dependence comes from the :1<->:2 dual-copy TI, NOT from
scaling the restraint." The run uses the v2.5.32 scaling method consistently, so
it is fine NOW. DANGER: if someone "fixes" the writer to match the v2.5.33 comment
(remove lambda scaling) WITHOUT switching the leg to real icfe=1 dual-copy TI,
the numerical restraint dG silently collapses toward 0 while the analytic term
still subtracts ~11.45 -> large silent error.
RECOMMENDED:
  1. Resolve the contradiction: either delete the stale v2.5.33 comment (scaling
     is the chosen method) OR finish the v2.5.33 migration (fixed k + icfe=1 TI).
  2. Add a consistency assertion at analysis time: the analytic correction's
     force constants MUST equal the lambda=1 window's actual &rst rk values
     (see boresch_consistency_assert.py).

## 12. Analysis path for the restraint leg (how to get the number)

The restraint leg is analyzed by the PACKAGE'S OWN analyzer (amber_md FEPAnalyzer,
fep.py): TI over per-window dV/dl (the DV/DL lines in each prod.out), with
alchemlyb MBAR/BAR as cross-checks used ONLY if they agree with TI within
2.0 kcal/mol (_ESTIMATOR_CONSISTENCY_KCAL). The leg total assembled in _leg()
(fep.py:1614-1632) = dG_kcal_mol (numerical TI) + dG_boresch_correction (-11.45).
The analyzer REQUIRES every window's prod.out; it returns None if any window is
missing (fep.py:1627) -> CANNOT run until lambda_0.400 finishes. The dt=0.001
rescue does NOT interfere: dV/dl is written to prod.out regardless of dt.
Do NOT hand-roll an estimator; run the package analyzer.

## 13. Artifacts produced this session (cc-file-api)

  * rescue_dt001_array.sh         -- parallel rescue, 5 windows (DONE, 5/5 PASS)
  * rescue_restraint_0400.sh      -- rescue the 6th window 0.400 (submitted)
  * audit_dt_all_legs.sh          -- read-only dt blast-radius audit
  * GENERATOR_FIX_restraint_dt.md -- the generator patch (restraint_dt_ps=0.001)
  * test_dt_fix_vdw_0100.sh       -- vdw dt A/B harness (SET ASIDE; no confirmed crash)
  * boresch_consistency_assert.py -- the latent-trap guard (sec 11a)
  * FEP_dt002_crash_investigation_HANDOFF.md -- this document
  * analyze_complex_restraint_FIXED.sh -- corrected analyzer submit (removes the
    stale v2.5.59 PYTHONPATH override; prints version provenance)
  * restraint_ti_crosscheck.py     -- numerical-vs-analytic restraint TI check (sec 14)

---

# ===== UPDATE 2026-06-23 (afternoon): RESTRAINT LEG ANALYZED & CLOSED =====

## 14. Restraint-leg analysis: result + TWO workflow gaps it exposed

The restraint leg is COMPLETE (all 20 windows show pmemd completion markers;
rescued 0.100-0.400 included). Attempting to analyze it surfaced two real,
distinct problems in the package -- documented here so they are not rediscovered.

### 14a. RESULT (what to report)
  dG_restraint = -11.45 kcal/mol  (the ANALYTIC Boresch standard-state
  correction; boresch_correction.txt = -11.447503; independently re-derived
  from on-disk lambda=1 &rst + geometry to within 0.01; package self-test
  reproduces Deng-Roux T4 11.62). This is the trustworthy restraint contribution.
  The SIMULATED restraint windows are NOT a usable TI source (see 14c).

### 14b. GAP 1 -- FEPAnalyzer CANNOT analyze the Option-A restraint leg
Running FEPAnalyzer(leg, lambdas).run() returned 0/20 windows, dG=None. Cause:
the analyzer's alchemlyb backend (extract_dHdl/extract_u_nk) requires the
Amber icfe=1 free-energy section (clambda / DV/DL). The restraint leg is
Option-A: icfe=0, lambda-scaled Boresch &rst, NO clambda, NO DV/DL line. Every
window failed with alchemlyb "No free energy section found ... clambda was None".
The non-alchemlyb fallback (_fallback_ti -> _parse_dvdl) ALSO found nothing,
because _parse_dvdl looks for a literal DV/DL token in the AVERAGES block that
icfe=0 never writes. So the analyzer is STRUCTURALLY unable to analyze this leg.
This is the v2.5.32/v2.5.33 schism (sec 11a) confirmed in practice: the
generator emits Option-A, the analyzer expects icfe=1 dual-copy TI. They were
never reconciled.
NOTE: this does NOT affect decharge/vdw/solvent legs -- those ARE icfe=1, so
FEPAnalyzer will work on them once they are complete from the clean run.

### 14c. GAP 2 -- the Option-A leg is not a sound NUMERICAL TI source
The only per-window alchemical signal is <RESTRAINT> in the AVERAGES block.
For U(lambda)=lambda*U_full: dU/dlambda = <RESTRAINT_lambda>/lambda, so
dG_attach = integral_0^1 <RESTRAINT>/lambda dlambda. That integrand has a
1/lambda SINGULARITY at lambda->0. On-disk <RESTRAINT> (decreases with lambda
because the fluctuation <(x-x0)^2> shrinks as k stiffens):
  | lambda | <RESTRAINT> | g=<R>/lambda |
  |---|---|---|
  | 0.100 | 34.35 | 343.5  <- singular driver |
  | 0.500 | 18.39 |  36.8 |
  | 1.000 | 10.18 |  10.2 |
restraint_ti_crosscheck.py integrates g three ways (3-point preview; full
20-window run refines but does not change the verdict):
  (a) naive trapezoid in lambda        : ~88 kcal/mol
  (b) sqrt(lambda) substitution        : ~63 kcal/mol
  (c) ln(lambda) + linear tail to 0    : ~87 kcal/mol
  analytic |correction|                : 11.45 kcal/mol
The ~40% spread among (a)/(b)/(c) IS the singularity artifact; the dominant
contribution lives in the UNSAMPLED lambda<0.1 region. None land near 11.45.
=> The numerical leg is unreliable BY CONSTRUCTION (linear restraint scaling
   endpoint singularity). Use the analytic term. The crosscheck script exists
   to DEMONSTRATE this, not to produce a reported number.
IMPORTANT for _leg(): _leg does total += rest['dG_kcal_mol'] + bcorr. If a
future analyzer ever returns a numerical dG_kcal_mol for this Option-A leg, it
would CORRUPT the total (double-count / add a singular number on top of the
analytic -11.45). The restraint leg's dG_kcal_mol should be 0 (or omitted) and
only the analytic correction applied.

### 14d. RECOMMENDED generator/analyzer fix (for the maintainer)
Pick ONE and make generator + analyzer agree:
  Option X (simplest, matches current physics): the restraint contribution is
    ANALYTIC ONLY (-11.45). Make the analyzer recognize an Option-A restraint
    leg (icfe=0 marker / absence of clambda) and return dG_kcal_mol=0 +
    dG_boresch_correction=<analytic>, skipping alchemlyb entirely. No singular
    TI. Cleanest.
  Option Y (if a simulated restraint dG is truly wanted): change the generator
    to a SOFT/nonlinear restraint schedule (e.g. lambda^k or a polynomial that
    flattens dU/dlambda near 0) AND densify lambda<0.1, AND write a matching
    icfe=1-style dHdl the analyzer can read. Much more work; rarely beats the
    analytic term for a Boresch restraint.
Either way: ADD the boresch_consistency_assert.py guard (sec 11a) so the
analytic term and the lambda=1 &rst constants can never silently diverge.

## 15. FINAL STATE / what remains for a full dG_bind

  COMPLEX restraint : DONE -> -11.45 kcal/mol (analytic), leg closed.
  COMPLEX decharge  : needs clean run; lambda_0.175 had died (dt=0.002). icfe=1
                      -> FEPAnalyzer works once complete.
  COMPLEX vdw       : needs clean run; killed-not-crashed (sec 10); diagnose the
                      ~0.600 region from REAL crash data after the clean run.
  SOLVENT decharge  : not yet run.
  SOLVENT vdw       : not yet run.
  Charge correction : auto (Rocklin) if ligand net-charged; self-skips if neutral.
  dG_bind = (complex: decharge + vdw + restraint[-11.45])
            - (solvent: decharge + vdw)  [+ charge corr if applicable]
NEXT ACTIONS:
  1. Apply GENERATOR_FIX_restraint_dt.md (restraint_dt_ps=0.001) -- and strongly
     consider extending dt=0.001 to ALL legs for THIS system (dt=0.002 proven on
     the stability edge here).
  2. Clean full run of all legs.
  3. Analyze decharge/vdw/solvent with the FIXED analyzer script (FEPAnalyzer
     works on those icfe=1 legs).
  4. Restraint contribution = -11.45 (do NOT re-run FEPAnalyzer on it; it cannot
     analyze Option-A -- see 14b).
  5. File GAP 1 + GAP 2 + sec 11a with the package maintainer.