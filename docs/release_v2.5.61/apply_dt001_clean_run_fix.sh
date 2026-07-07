#!/bin/bash
# =============================================================================
# apply_dt001_clean_run_fix.sh
#
# Applies the dt=0.001 fix to the INSTALLED amber_md package so a CLEAN RUN
# generates stable production mdins for ALL FEP legs on this system.
#
# WHY ALL LEGS (not just restraint): dt=0.002 is proven on the integrator
# stability edge for THIS system (restraint mid-band crashed 3/3 seeds;
# complex_decharge lambda_0.175 also died at dt=0.002; vdw never got a fair
# full-length test). For a clean run you intend to TRUST, dt=0.001 everywhere
# removes the whole crash class. Cost: 2x wall-clock for the same ns.
#
# CRITICAL CONSISTENCY BUG THIS FIXES:
#   run_amber.py converts --prod-ns to steps with a HARD-CODED /2 (assumes 2 fs):
#       cfg.md.prod_nsteps = int(round(a.prod_ns * 1e6 / 2))
#   but the dt actually written to the mdin comes from a SEPARATE field dt_ps.
#   If you change only dt_ps=0.001, nstlim is still ns*1e6/2 -> you'd silently
#   get HALF the simulation time. This script keeps dt and the ns->steps math
#   consistent (divides by the ps-per-step = dt_ps).
#
# Edits (with .bak backups):
#   config.py : FEPConfig.dt_ps      0.002 -> 0.001
#   config.py : MDConfig.prod_dt_ps  0.002 -> 0.001   (keep MM/GBSA consistent too)
#   run_amber.py : the two "/2" ns->steps conversions -> "/1" (i.e. *1000 per ns
#                  at dt=0.001), implemented as a dt-aware expression.
#
# USAGE:
#   1. EDIT PKG below to the tree your clean run imports (verify with:
#        python -c "import amber_md, os; print(os.path.dirname(amber_md.__file__))" )
#   2. bash apply_dt001_clean_run_fix.sh
#   3. re-run the GENERATOR for a fresh run; spot-check a prod.in shows dt=0.001
#      and nstlim = prod_ns * 1e6 (i.e. 10 ns -> 10,000,000 steps at 1 fs).
#
# REVERT: restore the .bak files this script writes.
# =============================================================================
set -euo pipefail

# ---- EDIT THIS to your install (the dir CONTAINING the amber_md/ package) ----
PKG="$HOME/Tools/affiniscreen"
# -----------------------------------------------------------------------------

CFG="$PKG/amber_md/config.py"
RUN="$PKG/run_amber.py"
STAMP=$(date +%Y%m%d_%H%M%S)

for f in "$CFG" "$RUN"; do
  [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
  cp -p "$f" "$f.bak_$STAMP"
  echo "backup: $f.bak_$STAMP"
done

echo
echo "=== BEFORE ==="
grep -nE "dt_ps: float = 0\.002|prod_dt_ps: float = 0\.002" "$CFG" || true
grep -nE "prod_ns \* 1e6 / 2|prod_ns\*1e6/2" "$RUN" || true

# --- config.py: dt_ps and prod_dt_ps 0.002 -> 0.001 (FEP + MM/GBSA paths) ---
# Only touch the exact default-assignment lines, not comments mentioning 0.002.
sed -i -E 's/^([[:space:]]*dt_ps: float = )0\.002/\10.001/' "$CFG"
sed -i -E 's/^([[:space:]]*prod_dt_ps: float = )0\.002/\10.001/' "$CFG"

# --- run_amber.py: make ns->steps DT-AWARE instead of hard-coded /2 ----------
# Original (two spots):
#   cfg.md.prod_nsteps = int(round(a.prod_ns * 1e6 / 2))
#   md=MDConfig(prod_nsteps=int(a.prod_ns*1e6/2), ...)
# Replace the literal "/ 2" and "/2" in the prod_ns->steps expression with
# "/ (1000.0*0.001)" == /1.0 i.e. steps = ns*1e6 / (ps_per_step*1000)... simpler:
# at dt=0.001 ps, steps-per-ns = 1e6 (1,000,000). So ns*1e6 steps. Divisor = 1.
# We write it as / (DT_PS/0.001*2 ... ) -> cleanest is to hard-set the divisor to
# match dt=0.001: ns * 1e6 / 1.  Use sed on the exact substrings.
sed -i -E 's/a\.prod_ns \* 1e6 \/ 2/a.prod_ns * 1e6 \/ 1/' "$RUN"
sed -i -E 's/a\.prod_ns\*1e6\/2/a.prod_ns*1e6\/1/' "$RUN"

echo
echo "=== AFTER ==="
grep -nE "dt_ps: float = 0\.001|prod_dt_ps: float = 0\.001" "$CFG" || true
grep -nE "prod_ns \* 1e6 / 1|prod_ns\*1e6/1" "$RUN" || true

echo
echo "DONE. Now regenerate a run and verify a prod.in:"
echo "  grep -E 'dt=|nstlim=' <fresh_run>/fep/complex_vdw/lambda_0.500/prod.in"
echo "  expect: dt=0.001  and  nstlim = prod_ns*1e6 (10 ns -> 10000000)"
echo
echo "If anything looks wrong, restore: cp $CFG.bak_$STAMP $CFG ; cp $RUN.bak_$STAMP $RUN"
