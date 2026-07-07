# v2.5.0 build final33: ABFE QC reader (Step 3)

## New: amber_md/abfe_qc.py
Post-run QC for an ABFE work-dir. Reads the REAL schema written by the
analyzer + cycle-closer (ABFE_RESULT.json + per-leg summary.json for
complex_decharge/complex_vdw/solvent_decharge/solvent_vdw).

Reports + verdict (GO / REVIEW): cycle closure & trusted flag; per-leg dG
breakdown; window completeness; MBAR adjacent overlap (warn <0.03);
estimator spread / convergence (warn >1.0 kcal/mol).

Usage (after the smoke test finishes):
  cd ~
  source $WF/activate_amber_md.sh
  python -m amber_md.abfe_qc ~/abfe_smoketest_<timestamp>
  python -m amber_md.abfe_qc ~/abfe_smoketest_<timestamp> --json

## Verified against synthetic results matching the exact schema:
  healthy -> GO (exit 0); broken -> REVIEW (exit 1) flags untrusted cycle,
  3/5 windows, weak overlap 0.005, spread 1.8.  All .py compile.
