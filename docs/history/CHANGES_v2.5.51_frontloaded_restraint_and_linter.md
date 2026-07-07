# v2.5.51 - front-loaded restraint ramp + run-script linter

Date: 2026-06-21

## restraint_lambdas_fine: 20-window quadratic front-load (opt-in)
Synthesized from FEP+ .msj (lig_12944901, converged; quadratic, ~38% <0.15), BAT.py attach_rest,
and the reference platform zero-shot. Grid: 0.0,0.004,0.016,0.036,0.064,0.1,0.144,0.196,0.256,0.324,0.4,0.5,0.6,0.7,
0.8,0.875,0.925,0.95,0.975,1.0. Default restraint_lambdas unchanged.

## tools/lint_run_script.py
R1 errexit+pipefail; R2 unguarded bare parses; R3 escalation markers; R4 array ranges.
