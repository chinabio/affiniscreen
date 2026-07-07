# amber_md v2.5.68

## Restraint MD leg OFF by default
- The explicit ~2 ns Boresch restraint MD leg is no longer run by default.
- Its sampled free energy is ~0 (removing a well-centered restraint on the
  bound, fully-interacting complex). The analytic Boresch standard-state
  correction is folded onto complex_vdw instead (Boresch 2003 / FEP+ / OpenFE).
- INVARIANT (CI-guarded): the Boresch POTENTIAL is still held ON during
  complex_decharge + complex_vdw, so the decoupled ligand cannot drift.
- `--restraint-leg` restores the explicit MD leg (BAT.py-style validation).
- analyze_campaign.py reads the analytic term from complex_vdw when no
  restraint leg is present; ABFE result is identical to the with-leg cycle
  when the restraint-apply FE is ~0.
