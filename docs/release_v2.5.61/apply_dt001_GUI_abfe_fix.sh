#!/bin/bash
# =============================================================================
# apply_dt001_GUI_abfe_fix.sh
#
# Applies the dt=0.001 fix for runs launched through the GUI (Streamlit
# Setup_and_Launch page), method = ABFE, engine = Amber -- i.e. the user's
# actual launch path. Verified against amber_md v2.5.60.
#
# WHY THE EARLIER run_amber.py PATCH WAS WRONG FOR YOU:
#   The GUI does NOT use run_amber.py's argparse defaults. For ABFE/Amber it
#   builds the command in gui/pages/0_Setup_and_Launch.py and:
#     * NEVER passes --dt  -> dt comes from config.py FEPConfig.dt_ps (=0.002)
#     * computes --nstlim-prod = complex_ns * 1e6 / 2  (HARD-CODED 2 fs)
#   So dt and nstlim are set in TWO disconnected places, both assuming 2 fs.
#   To move to 1 fs WITHOUT silently halving the simulated time, BOTH must change.
#
# WHAT THIS PATCHES (with .bak backups):
#   config.py:130  dt_ps: float = 0.002        -> 0.001   (the dt ABFE actually uses)
#   config.py:44   prod_dt_ps: float = 0.002   -> 0.001   (MD/heat path consistency)
#   0_Setup_and_Launch.py  ABFE branch (~708-713): _cns*1e6/2 -> _cns*1e6/1
#                                                   _ens*1e6/2 -> _ens*1e6/1
#   0_Setup_and_Launch.py  RBFE branch (~549-554): same (only matters if you run RBFE)
#   fep_common.py:230  dt widget default 0.002  -> 0.001   (defensive; relative-FEP UI)
#
# RESULT: a GUI ABFE/Amber launch with complex_ns=10 now generates per-window
#   prod.in with dt=0.001 and nstlim=10,000,000 (= 10 ns at 1 fs), for ALL legs
#   (restraint + decharge + vdw + solvent). dt=0.002 is proven on the stability
#   edge for this system, so all-legs dt=0.001 removes the whole crash class.
#   Cost: 2x wall-clock for the same ns.
#
# USAGE:
#   1. verify PKG below matches:
#        python -c "import amber_md, os; print(os.path.dirname(amber_md.__file__))"
#      (should print  <PKG>/amber_md )
#   2. bash apply_dt001_GUI_abfe_fix.sh
#   3. RESTART the Streamlit GUI (so it re-imports the patched modules).
#   4. Launch a fresh ABFE/Amber run; before committing GPU, spot-check ONE prod.in:
#        grep -E 'dt=|nstlim=' <run>/lig_*/fep/complex_vdw/lambda_0.500/prod.in
#        expect:  dt=0.001   and   nstlim=10000000  (for complex_ns=10)
#
# REVERT: restore the .bak_<stamp> files this script writes.
# =============================================================================
set -euo pipefail

# ---- EDIT if needed: dir CONTAINING the amber_md/ package -------------------
PKG="$HOME/Tools/affiniscreen"
# -----------------------------------------------------------------------------

CFG="$PKG/amber_md/config.py"
SETUP="$PKG/amber_md/gui/pages/0_Setup_and_Launch.py"
FEPC="$PKG/amber_md/gui/fep_common.py"
STAMP=$(date +%Y%m%d_%H%M%S)

for f in "$CFG" "$SETUP" "$FEPC"; do
  [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }
  cp -p "$f" "$f.bak_$STAMP"; echo "backup: $f.bak_$STAMP"
done

echo; echo "=== BEFORE ==="
grep -nE "dt_ps: float = 0\.002|prod_dt_ps: float = 0\.002" "$CFG" || true
grep -nE "_cns \* 1e6 / 2|_ens \* 1e6 / 2" "$SETUP" || true
grep -nE 'number_input\("dt \(ps\)", 0\.0005, 0\.004, 0\.002' "$FEPC" || true

# --- config.py: dt defaults 0.002 -> 0.001 (exact default-assignment lines) ---
sed -i -E 's/^([[:space:]]*dt_ps: float = )0\.002/\10.001/'      "$CFG"
sed -i -E 's/^([[:space:]]*prod_dt_ps: float = )0\.002/\10.001/' "$CFG"

# --- GUI nstlim computation: /2 -> /1  (ONLY the _cns/_ens ->steps lines) -----
# These are the complex_ns / equil_ns -> nstlim conversions in BOTH the RBFE
# (549/554) and ABFE (708/713) branches. Matching the exact substrings avoids
# touching any other 1e6 usage.
sed -i -E 's/(_cns \* 1e6 \/ )2\)/\11)/g' "$SETUP"
sed -i -E 's/(_ens \* 1e6 \/ )2\)/\11)/g' "$SETUP"

# --- fep_common.py: dt widget default 0.002 -> 0.001 (defensive) --------------
sed -i -E 's/(number_input\("dt \(ps\)", 0\.0005, 0\.004, )0\.002/\10.001/' "$FEPC"

echo; echo "=== AFTER ==="
grep -nE "dt_ps: float = 0\.001|prod_dt_ps: float = 0\.001" "$CFG" || true
grep -nE "_cns \* 1e6 / 1|_ens \* 1e6 / 1" "$SETUP" || true
grep -nE 'number_input\("dt \(ps\)", 0\.0005, 0\.004, 0\.001' "$FEPC" || true

echo
echo "DONE. RESTART the GUI, launch a fresh ABFE/Amber run, then verify ONE prod.in:"
echo "  grep -E 'dt=|nstlim=' <run>/lig_*/fep/complex_vdw/lambda_0.500/prod.in"
echo "  expect: dt=0.001  and  nstlim=10000000  (complex_ns=10)"
echo
echo "REVERT: cp <file>.bak_$STAMP <file>  for each of:"
echo "  $CFG"
echo "  $SETUP"
echo "  $FEPC"
