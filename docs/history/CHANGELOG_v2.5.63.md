# amber_md v2.5.63  (builds on v2.5.62)

## Fixed
- Option-A restraint leg now analyzable. complex_restraint is plain NPT MD
  (prod.in has no clambda/icfe); alchemlyb/MBAR reported "no free energy
  section" and analyze exited 2 (dG=None). New fep.analyze_restraint_leg_optionA()
  detects this mode, scores the leg from the analytic Boresch correction, writes a
  COMPLETE summary.json, exits 0. Real TI restraint legs still use FEPAnalyzer/MBAR.

## Changed
- True 10 ns default: fep_driver --nstlim-prod 5_000_000 -> 10_000_000 so the CLI
  matches config (10 ns @1fs) instead of overriding it with 5 ns. (The restraint
  leg does not need 10 ns -- ~1 ns suffices for an analytic correction -- this is
  for cross-leg consistency and to remove the override footgun.)

## Carried from v2.5.62
- dt 0.002->0.001 everywhere; step counts x2; _dt_guard wired into _build_configs;
  tools/check_dt_regression.py (now also asserts 10 ns default + restraint analyzer).
