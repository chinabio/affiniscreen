# v2.5.44 - pipefail parse hotfix (silent prod abort)

Date: 2026-06-21

Symptom (v2.5.43 run lig_12944901): min/heat/dens/eq PASS incl. boresch gate, but no prod.out / prod.console / recovery.log on any window; analyzer 0/16. The job exited right after 'equilibrated geometry OK', never printing 'prod attempt 1/10'.

Root cause: script runs under 'set -euo pipefail'. In run_prod_with_restart the lines
  _CL=$(grep -oE 'clambda=...' prod.in | head -1 | cut -d= -f2)   (and timask1/2, ifsc, icfe)
are BARE command-substitution assignments. The Option A prod.in contains none of those keywords, so grep exits 1; pipefail propagates rc=1; and because the assignment is not prefixed with 'local' or guarded, set -e aborts the window before the first pmemd prod attempt -- with no recovery log.

Fix: append '|| true' to all seven prod.in parses so a no-match yields an empty string instead of aborting. The downstream ': "${_CL:=0.0}"' / ': "${_IFSC:=0}"' defaults already handle empties, and the icfe-gated clambda guard (v2.5.41) still catches a genuinely malformed TI prod.in.
