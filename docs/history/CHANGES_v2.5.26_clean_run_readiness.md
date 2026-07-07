# v2.5.26 -- clean-run readiness + version track bump

- Heat-stage stability gate in both run paths (peak T > heat_temp_max_K -> exit 70 before dens).
- tools/verify_stage_order.sh GPU-free pre-submit dry-checker.
- README_CLEAN_RUN.md clean-run checklist.
- VERSION TRACK: __version__ 2.5.23 -> 2.5.26 (single source of truth in
  amber_md/__init__.py); __build__ -> 'clean-run-ready'; VERSION file -> 2.5.26;
  run_amber.py __version__ 2.5.22 -> 2.5.26 (was stale); changelog entries added
  for 2.5.24, 2.5.25, 2.5.26. The .lsf banner (amber_md.version.lsf_banner) now
  stamps v2.5.26 into every generated job script.
