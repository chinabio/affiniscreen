# v2.5.45 - full pipefail audit of generated scripts

Date: 2026-06-21

Follow-up to v2.5.44. Systematically scanned every shell assignment of the form VAR=$(pipeline) in
fep.py-generated scripts that run under 'set -euo pipefail'. Findings:
  - 7 prod.in parses in run_prod_with_restart: FIXED in v2.5.44 (|| true).
  - eq/heat peak-TEMP scans (maxT/maxTh): already safe -- wrapped in 'set +o pipefail; set +e; ...; true'.
  - STAGE FAILED grep | tail handlers: already guarded with '|| true'.
  - HREMD per-window crash scan 'hit=$(grep ... | tail -n 3)': safe-by-accident (tail last in pipe);
    hardened in v2.5.45 with explicit '|| true'.
No other unguarded bare command-substitution assignments remain.
