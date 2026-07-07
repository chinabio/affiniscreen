#!/usr/bin/env bash
# run_screen_analysis.sh -- walk every lig_* in a screen directory,
# then build the master ranking page.
set -euo pipefail
SCREEN=${1:?"usage: $0 <screen_dir>"}
KIT=$(dirname "$(readlink -f "$0")")

shopt -s nullglob
for LIG in "$SCREEN"/lig_*; do
  [[ -d "$LIG" ]] || continue
  echo "==> $(basename "$LIG")"
  "$KIT/run_analysis.sh" "$LIG" || echo "  !! failed, continuing"
done

# Master ranking page
python3 "$KIT/build_screen_summary.py" "$SCREEN"

echo
echo "Done."
echo "  Per-ligand reports: <ligand>/analysis/COMBINED_REPORT.html"
echo "  Screen summary:     $SCREEN/screen_summary.html"
