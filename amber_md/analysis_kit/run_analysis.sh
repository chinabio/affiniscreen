#!/usr/bin/env bash
# run_analysis.sh  (v5) -- per-ligand post-processing, path-configurable.
#
# Usage:
#   ./run_analysis.sh [--dry-run] <ligand_workdir>
#
# Configurable via env-vars (paths RELATIVE to <ligand_workdir>):
#   TOP_REL      default: auto-detect, preferring solvated/complex_solv > complex
#   TRAJ_REL     default: auto-detect, preferring prod*.nc
#   LIG_RESNAME  default: LIG
set -euo pipefail

DRY=0
if [[ "${1:-}" == "--dry-run" ]]; then DRY=1; shift; fi

LIG=${1:?"usage: $0 [--dry-run] <ligand_workdir>"}
LIG=$(readlink -f "$LIG")
KIT=$(dirname "$(readlink -f "$0")")

# ----------------------------------------------------------------------
# Topology auto-detect.
# Trajectory was integrated against the FULL SOLVATED system, so we MUST
# pick the solvated topology. Order of preference:
#   1) solvated.prmtop / solvated.parm7
#   2) complex_solv.{prmtop,parm7}    (other naming conventions)
#   3) system.prmtop  (last well-known name)
#   4) anything that ISN'T receptor/ligand/complex-without-solvent
# ----------------------------------------------------------------------
detect_top () {
  local base="$1" n cand
  for n in solvated.prmtop solvated.parm7 \
           complex_solv.prmtop complex_solv.parm7 \
           system.prmtop system.parm7; do
    [[ -f "$base/topo/$n" ]] && { echo "topo/$n"; return; }
    [[ -f "$base/$n"      ]] && { echo "$n";      return; }
  done
  # Heuristic fallback: any prmtop NOT named like the small/intermediate ones
  cand=$(find "$base" -maxdepth 3 \( -name '*.prmtop' -o -name '*.parm7' \) \
          | grep -Ev '/(ligand|lig|small|receptor|apo|complex)\.(prmtop|parm7)$' \
          | head -n1 || true)
  [[ -n "$cand" ]] && echo "${cand#$base/}" && return
  # Last resort
  cand=$(find "$base" -maxdepth 3 \( -name '*.prmtop' -o -name '*.parm7' \) | head -n1 || true)
  [[ -n "$cand" ]] && echo "${cand#$base/}"
}

detect_traj () {
  local base="$1" cand
  cand=$(find "$base" -maxdepth 5 -type f \( -name 'prod*.nc' -o -name 'production*.nc' \
           -o -name 'prod*.dcd' \) 2>/dev/null | sort | head -n1 || true)
  [[ -n "$cand" ]] && echo "${cand#$base/}" && return
  cand=$(find "$base" -maxdepth 5 -type f \( -name '*.nc' -o -name '*.dcd' \) 2>/dev/null \
          | grep -Ev '/(min|heat|equil|nvt|npt|cool|warm)[^/]*\.(nc|dcd)$' \
          | sort | head -n1 || true)
  [[ -n "$cand" ]] && echo "${cand#$base/}" && return
  cand=$(find "$base" -maxdepth 5 -type f \( -name '*.nc' -o -name '*.dcd' \) 2>/dev/null | head -n1 || true)
  [[ -n "$cand" ]] && echo "${cand#$base/}"
}

TOP_REL=${TOP_REL:-}
TRAJ_REL=${TRAJ_REL:-}
LIG_RESNAME=${LIG_RESNAME:-LIG}

if [[ -z "$TOP_REL" || ! -f "$LIG/$TOP_REL" ]]; then
  [[ -n "$TOP_REL" ]] && echo "[warn] TOP_REL='$TOP_REL' not found, auto-detecting..." >&2
  TOP_REL=$(detect_top "$LIG" || true)
fi
[[ -z "$TOP_REL" || ! -f "$LIG/$TOP_REL" ]] && {
  echo "ERROR: no topology found under $LIG (set TOP_REL=...)" >&2
  find "$LIG" -maxdepth 4 \( -name '*.prmtop' -o -name '*.parm7' \) >&2
  exit 1
}

if [[ -z "$TRAJ_REL" || ! -f "$LIG/$TRAJ_REL" ]]; then
  [[ -n "$TRAJ_REL" ]] && echo "[warn] TRAJ_REL='$TRAJ_REL' not found, auto-detecting..." >&2
  TRAJ_REL=$(detect_traj "$LIG" || true)
fi
[[ -z "$TRAJ_REL" || ! -f "$LIG/$TRAJ_REL" ]] && {
  echo "ERROR: no trajectory found under $LIG (set TRAJ_REL=...)" >&2
  find "$LIG" -maxdepth 5 -type f \( -name '*.nc' -o -name '*.dcd' \) >&2 || true
  exit 1
}

# Warn if we somehow ended up with a tiny topology (likely ligand- or apo-only)
top_bytes=$(stat -c%s "$LIG/$TOP_REL" 2>/dev/null || echo 0)
if [[ "$top_bytes" -lt 1000000 ]]; then
  echo "[warn] TOP_REL='$TOP_REL' is small (${top_bytes} B); analysis needs the SOLVATED topology." >&2
fi

echo "[run_analysis] LIG=$LIG"
echo "[run_analysis]   TOP_REL  = $TOP_REL"
echo "[run_analysis]   TRAJ_REL = $TRAJ_REL"
echo "[run_analysis]   LIG_RES  = $LIG_RESNAME"

if [[ "$DRY" -eq 1 ]]; then
  echo "[run_analysis] --dry-run requested; not running cpptraj."
  exit 0
fi

mkdir -p "$LIG/analysis"
cd       "$LIG/analysis"

TOP_FROM_ANALYSIS="../$TOP_REL"
TRAJ_FROM_ANALYSIS="../$TRAJ_REL"

render () {
  local src="$1" dst="$2"
  sed -e "s|{{TOP_REL}}|$TOP_FROM_ANALYSIS|g" \
      -e "s|{{TRAJ_REL}}|$TRAJ_FROM_ANALYSIS|g" \
      -e "s|{{LIG}}|$LIG_RESNAME|g" \
      "$src" > "$dst"
}

render "$KIT/cpptraj_rmsd_rmsf.in.template" cpptraj_rmsd_rmsf.in
render "$KIT/load_pymol.pml.template"       load_pymol.pml
render "$KIT/load_vmd.tcl.template"         load_vmd.tcl

cpptraj -i cpptraj_rmsd_rmsf.in > cpptraj.log 2>&1 || {
  echo "ERROR: cpptraj failed. See $LIG/analysis/cpptraj.log" >&2
  tail -30 cpptraj.log >&2
  exit 1
}

python3 "$KIT/combine_report.py" "$LIG"
echo "[run_analysis] done -> $LIG/analysis/COMBINED_REPORT.html"
