# CHANGES - amber_md ABFE correctness release

## v2.5.19  (2026-06-16)

Cleanup + equilibration hardening on top of the v2.5.18 correctness wiring.

### Removed
* amber_md/fep_mdin_v2516.py -- dead module (zero references anywhere). Its
  role (validated GTI block + lambda schedules) lives in config.py/fep.py.

### Equilibration: positional restraints through dens + eq (patch sec.4)
* fep.py: new _posres_block() emits ntr/restraint_wt/restraintmask for the
  density and equilibration stages, holding solute heavy atoms while the box
  and soft-core ligand relax -- the FEP-SPell staged-equilibration behaviour.
  Wired into _dens_in (stiff, min_restraint_wt) and _eq_in (gentle,
  equil_restraint_wt). Uses the EXISTING min->dens->eq->prod runner, so the
  hardened restart / box-drift / temperature-gate shell logic is untouched.
* config.py: FEPConfig.posres_mask (default "" -> unchanged behaviour; set it
  to a heavy-atom mask to enable, e.g. "!:WAT,Na+,Cl-,K+ & !@H=").

### Intentionally KEPT (removal would break imports / lose the reference)
* boresch_restraint.py  -- canonical FEP-SPell reference the production writer
  is validated against (6/6 DOF self-consistency test).
* boresch_autoselect.py -- imported by fep.py (production).
* abfe_integration.py / equilibration_fepspell.py -- public helpers used by
  the smoke test; not on the production path but harmless to keep.

### Equilibration engine NOT rerouted (deliberate)
* equilibration_fepspell.write_equilibration emits a variable multi-step ladder
  driven by RUN_ORDER.txt. The production runner hard-codes min->dens->eq->prod
  with failure recovery; rerouting would require rewriting that shell logic for
  no thermodynamic benefit (dens/eq are already gentle Berendsen/dt=0.001).
  Closed the only real gap (carried positional restraints) in place instead.

---

# CHANGES - amber_md ABFE correctness release

## v2.5.18  (2026-06-16)

Wires the FEP-SPell-ABFE-faithful restraint, analytic correction, and charge
correction into the PRODUCTION runtime path (fep.py / fep_driver.py). The
v2.5.17 release had ported these modules but left them orphaned -- only the
smoke test imported them, while the code that actually generates and submits
windows still used the legacy, self-inconsistent implementations.

### Restraint self-consistency (THE correctness fix)
* fep.py._write_boresch_RST now emits the six Boresch &rst records in the exact
  FEP-SPell-ABFE atom ordering:
      r=L1-P1, alpha=P1-L1-L2, theta=P2-P1-L1,
      gamma=P1-L1-L2-L3, beta=P2-P1-L1-L2, phi=P3-P2-P1-L1
  Previously the writer restrained a DIFFERENT set of angles/dihedrals than the
  workflow measured and corrected (only 1/6 DOF matched), biasing every ABFE.
  Now the SIMULATED, MEASURED, and CORRECTED coordinates are identical (6/6).
* boresch.py: added canonical_dofs_from_legacy() mapping the legacy
  select_boresch_atoms dict (aA/bA/cA + A/B/C, thA0/thB0/phA0/phB0/phC0) to the
  canonical FEP-SPell six DOF (verified numerically: alpha=thB0, theta=thA0,
  gamma=phC0, beta=phB0, phi=phA0).

### Analytic standard-state correction (Deng & Roux Eq.38/40)
* boresch.py: added restraint_correction_dengroux() and
  boresch_correction_dengroux() (verbatim FEP-SPell analytic.py, rk->2*rk).
  Reproduces the published 11.62 kcal/mol reference (computed 11.618).
* fep_driver.py: the ABFE correction now uses boresch_correction_dengroux()
  evaluated on the SAME geometry/force constants written to boresch.RST, for
  both auto-selected atoms and a boresch.json loaded from disk. The legacy
  closed-form correction (evaluated on an inconsistent geometry) is no longer
  used for ABFE.

### Charge correction (Rocklin PB finite-size) wired in
* fep.py.build_cycle_closer_lsf now AUTOMATICALLY computes the
  complex-minus-solvent finite-size correction for net-charged ligands from the
  final frame of each decoupled-end-state leg, writes charge_correction.json,
  and folds it into dG_bind. Self-skips (0.0) for neutral ligands or if any
  input/dependency is unavailable -- never aborts the cycle-closer.
* config.py: FEPConfig gained solvent_mask and water_model (consumed above).

### Carried forward from v2.5.16/2.5.17 (unchanged, verified retained)
* Complete GTI soft-core block; scalpha 0.2 / scbeta 50; gti_chg_keep 1.
* 44-window vdW grid (dense 0.575-0.80); decharge + 16-window restraint legs.
* Boresch rk->2*rk factor fix; setup-time + in-job geometry gates.

### Validation
* All package .py files parse.
* Deng-Roux dG_release = 11.618 kcal/mol (FEP-SPell ref 11.62).
* End-to-end self-consistency: fep.py-simulated restraint == FEP-SPell-ported
  gen_restraint_str == analytic correction geometry (6/6 DOF).
