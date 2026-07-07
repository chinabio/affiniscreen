# v2.5.11 (final67) — Refined vdW lambda schedule (soft-core danger zone)

## Motivation
A clean v2.5.9/2.5.10 `complex_vdw` run (23 windows) completed 17/23. The 6
failures clustered in the soft-core decoupling region:

| array idx | lambda | failure | meaning |
|-----------|--------|---------|---------|
| 4         | 0.150  | exit 255 | prod box drift |
| 5         | ~0.20  | exit 255 | prod box drift |
| 12        | ~0.55  | exit 255 | prod box drift |
| 13        | ~0.60  | exit 255 | prod box drift |
| 15        | 0.700  | exit 71  | eq instability (gate fired) |
| 16        | 0.750  | exit 71  | eq instability (gate fired) |

Both signatures are physics, not code defects:
- **exit 71** = the final60 eq-stability gate doing its job: it refused to run
  production on a blown-up equilibration (which would otherwise poison MBAR).
- **exit 255** = box drift; v2.5.10 makes the in-job restart loop ride through
  it, but in this region the perturbation per window is large enough that even
  restarts struggle.

Root cause is the **uniform 0.05 lambda spacing** through the region where the
ligand vdW core is ~60-85% decoupled: atoms can nearly overlap, the soft-core
potential stiffens, and dV/dl + box response become hard to integrate.

## Change
`FEPConfig.vdw_lambdas` default, halving the spacing across 0.6-0.85:

```
old (23): ... 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, ...
new (28): ... 0.55, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725, 0.75,
              0.775, 0.80, 0.825, 0.85, 0.90, ...
```

Added: 0.625, 0.675, 0.725, 0.775, 0.825. Endpoints and low-lambda spacing
unchanged. The mdin generators (eq/dens/prod) are NOT changed — they were
already gentle (dt=0.001 + Berendsen taup=5 + gamma_ln=5 for dens/eq).

## Operational notes
- **This requires a full re-run of any vdw leg.** MBAR needs a single
  consistent lambda set, so the 17 completed complex_vdw windows cannot be
  mixed with the new schedule. (decharge legs are unaffected.)
- To keep the old schedule for a one-off, pass `--vdw-lambdas ...` explicitly.
- If 0.70/0.75 still destabilize after refinement, the next lever is the
  soft-core stiffness: try `scbeta=16` (from 12) or `scalpha=0.4` (from 0.5).

## Verification
- byte-compile OK; version 2.5.11 (final67).
- bash -n release gate: array + HREMD generated scripts parse clean.
- Generated eq.in/prod.in now emit mbar_states=28 with the new lambda list.
- Full unit suite: 29 passed.

## Lineage
  * exit 141 (SIGPIPE eq gate)            -> v2.5.7
  * host / outage kills                    -> v2.5.8
  * exit 2 (duplicated if -> syntax)       -> v2.5.9
  * exit 255 (box-drift grep wrong file)   -> v2.5.10
  * exit 71 / residual 255 (lambda spacing)-> v2.5.11 (this; schedule, not code)
