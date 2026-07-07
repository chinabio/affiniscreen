# v2.5.0 build final41: vdW endpoint schedule (lambda->1 singularity)

## Why
The smoke test blew up at the bare vdW endpoint lambda=1.000:
    Energy at 0.0000 = ****************   (end-state softcore singularity)
User asked to "cap vdW at 0.95". A HARD cap (dropping 1.0) avoids the blow-up
but BIASES dG_decouple by the missing 0.95->1.0 contribution -- the ligand is
never fully decoupled. That is scientifically wrong for ABFE.

## What (config.py default vdw_lambdas)
Keep FULL decoupling but add fine spacing near the endpoint so each near-1.0
dV/dl step is small and numerically stable:
  before: ... 0.85, 0.9, 0.95, 1.0                       (21 windows)
  after : ... 0.85, 0.9, 0.95, 0.975, 0.99, 1.0          (23 windows)
  endpoint gaps now 0.05 / 0.025 / 0.015 / 0.01 (was a single 0.05 jump to 1.0)

This is the standard production cure for endpoint singularities (dense lambda
near the decoupled state) and pairs with final39's dV/dl guard + final40's
per-lambda decorrelation / honest MBAR handling.

## Hard 0.95 cap (if ever wanted)
Pass an explicit --vdw-lambdas 0.0 ... 0.95 and ACCEPT the truncated-tail bias.
Not the default.

## Verified
  new schedule: 23 windows, strictly monotonic, ends at 1.0, max gap 0.05,
  endpoint spacing <= 0.01. config.py compiles.
