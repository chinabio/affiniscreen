# amber_md v2.5.64  (builds on v2.5.63)

## Changed -- restraint leg gets its own (shorter) production length
The Option-A complex_restraint leg is plain NPT equilibration whose free-energy
contribution is the ANALYTIC Boresch correction; it is NOT sampling-limited, so
running it as long as the TI legs (10 ns) wastes time without improving the
result. New dedicated default:

    config.FEPConfig.restraint_nstlim_prod = 2_000_000   # 2 ns @ 1 fs

- _prod_in_restraint() now emits this value (falls back to nstlim_prod if 0).
- New CLI flag: --restraint-nstlim-prod (default 2_000_000), threaded into FEPConfig.
- TI legs (decharge/vdw/solvent) are unchanged at nstlim_prod (10 ns @ 1 fs).

Rationale: a Boresch-restrained complex relaxes in well under 1 ns; 2 ns is a
safety margin. This keeps the correct asymmetry: TI legs need long sampling,
the restraint leg does not.

## Carried from v2.5.63 / v2.5.62
- Option-A restraint analyzer (analytic Boresch -> complete summary, exit 0).
- True 10 ns TI default; dt=0.001 everywhere; _dt_guard; CI check.
