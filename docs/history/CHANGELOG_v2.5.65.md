# amber_md v2.5.65  (builds on v2.5.64)

## Added -- tools/analyze_leg.py (stand-alone leg re-analyzer)
Re-analyze an already-finished leg directory without LSF or re-simulation:

    python tools/analyze_leg.py <leg_dir>
    python tools/analyze_leg.py <leg_dir> --lambdas 0.0,0.15,...,1.0 --temp 298.0

- Runs the SAME logic as the generated analyze_<leg>.lsf: Option-A restraint
  leg (no clambda) -> analytic Boresch value; TI leg -> FEPAnalyzer/MBAR.
- Auto-discovers the lambda schedule from the leg's analyze_*.lsf, else from
  the lambda_* sub-directories (override with --lambdas).
- Writes summary.json and uses the pipeline's exit codes:
  0=OK, 1=incomplete, 2=no dG, 3=exception, 4=usage/discovery error.

## Carried
- v2.5.64 restraint_nstlim_prod (2 ns restraint vs 10 ns TI).
- v2.5.63 Option-A analytic analyzer; v2.5.62 dt=1fs fix + guard + CI check.
