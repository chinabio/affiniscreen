#!/usr/bin/env bash
# verify_stage_order.sh (v2.5.26) -- dry-check a generated leg BEFORE submitting.
set -u
leg="${1:-.}"; fail=0
note(){ printf '  %-7s %s\n' "$1" "$2"; }
echo "verify_stage_order: scanning $leg"
shopt -s nullglob
wins=("$leg"/lambda_*)
[ ${#wins[@]} -eq 0 ] && { echo "ERROR: no lambda_* under $leg" >&2; exit 2; }
for wd in "${wins[@]}"; do
  lam=$(basename "$wd"); miss=""
  for f in min.in dens.in eq.in prod.in; do [ -f "$wd/$f" ] || miss="$miss $f"; done
  [ -n "$miss" ] && { note FAIL "$lam missing:$miss"; fail=1; } || note ok "$lam mdin present"
  if [ -f "$wd/heat.in" ]; then
    grep -qE 'ntp *= *0' "$wd/heat.in" || { note FAIL "$lam heat.in not NVT"; fail=1; }
    grep -qE 'nmropt *= *1' "$wd/heat.in" || { note FAIL "$lam heat.in no TEMP_0 ramp"; fail=1; }
  fi
done
rs=$(ls "$leg"/run*.sh "$leg"/*.lsf "$leg"/submit*.sh 2>/dev/null | head -n1)
if [ -n "${rs:-}" ] && [ -f "$rs" ]; then
  echo "checking run-script order in: $rs"
  if grep -q 'heat.in' "$rs"; then
    lmin=$(grep -n -- '-i min.in'  "$rs" | head -n1 | cut -d: -f1)
    lheat=$(grep -n -- '-i heat.in' "$rs" | head -n1 | cut -d: -f1)
    ldens=$(grep -n -- '-i dens.in' "$rs" | head -n1 | cut -d: -f1)
    if [ -n "$lmin" ] && [ -n "$lheat" ] && [ -n "$ldens" ] && [ "$lmin" -lt "$lheat" ] && [ "$lheat" -lt "$ldens" ]; then
      note ok "run-script order: min -> heat -> dens"
    else note FAIL "order wrong (min=$lmin heat=$lheat dens=$ldens)"; fail=1; fi
    grep -q 'stage=heat INSTABILITY' "$rs" && note ok "heat gate present" || { note FAIL "heat gate missing"; fail=1; }
    grep -q 'DENS_C' "$rs" && note ok "dens consumes \$DENS_C" || { note FAIL "dens not wired to heat.rst"; fail=1; }
  else note ok "no heat.in (do_heat off) -- min -> dens"; fi
else note WARN "no run script under $leg"; fi
echo
[ "$fail" -eq 0 ] && echo "PASS: $leg is ready to submit" || { echo "FAIL: fix issues above" >&2; }
exit "$fail"
