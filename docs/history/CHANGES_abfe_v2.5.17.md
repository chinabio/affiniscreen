# CHANGES — amber_md ABFE correctness release

## v2.5.17  (2026-06-16)

Makes a clean ABFEP run fully automatic and self-consistent. Audited the ENTIRE
FEP-SPell-ABFE codebase and ported every piece needed between "topology built"
and "windows submitted".

### Correctness (restraint self-consistency) — IMPORTANT
* Boresch restraint atom ORDERING corrected. v2.5.16 restrained a different set
  of angles/dihedrals than it measured and analytically corrected. Now matches
  FEP-SPell exactly: r=L1-P1, alpha=P1-L1-L2, theta=P2-P1-L1,
  gamma=P1-L1-L2-L3, beta=P2-P1-L1-L2, phi=P3-P2-P1-L1. The simulated potential,
  the cpptraj measurement, and the analytic correction now use the SAME six
  internal coordinates.
* NEW six-DOF measurement (amber_md.boresch_restraint.measure_six_dofs): runs
  the FEP-SPell cpptraj template on the equilibrated trajectory, circular-means
  the dihedrals, feeds BOTH restraints.inp and the analytic correction.

### Automation (no manual atom picking required) — NEW
* amber_md.boresch_autoselect: automatic six-atom selection. RMSF-aware, Calpha
  receptor anchors in a 5-11 A shell, scored combinatorial search driving
  r/alpha/theta to well-conditioned values; rejects collinear/too-close picks.
  - ON BY DEFAULT; manual lig_/rec_restraint_atoms override it verbatim.
  - Always writes boresch_atom_selection.dat (QC report).
  - Raises with an actionable message if no good set exists -- never the silent
    'pass' restraint upstream FEP-SPell writes.

### Protocol (window stability) — NEW
* amber_md.equilibration_fepspell: staged equilibration ported verbatim --
  2-step minimization (restrained -> free) + 5->100->200->298.15 K heat ladder
  (3 heat + 3 press cycles, TEMP_0 ramp) + free relax.

### Integration
* amber_md.abfe_integration: resolve_restraint_atoms() (manual-or-auto),
  prepare_complex_boresch() (measure + write + correct),
  write_staged_equilibration().

### Carried forward from v2.5.16 (unchanged)
* Complete GTI soft-core block; scalpha 0.2 / scbeta 50; gti_chg_keep 1.
* 44-window vdW grid (dense 0.575-0.80); 15-window decharge; 16-window
  restraint-removal leg; Boresch rk->2*rk fix (11.629 vs ref 11.62); Rocklin PB
  charge correction; cycle-closer folds in MD restraint + analytic + charge.

### Validation
* abfe_smoke_test.py: 8/8 checks (added [8/8] integration) -- PASSED.
* boresch_autoselect: 8/8 random geometries yield a well-conditioned set.

### v2.5.17 addendum — post-equilibration restraint gate
* abfe_integration.verify_or_reselect_boresch(): after equilibration, validates
  the six Boresch atoms on the EQUILIBRATED structure. If invalid: AUTO picks are
  re-selected on the equilibrated coords (then progressively relaxed thresholds);
  MANUAL picks fail fast (never silently overridden). Raises with a clear message
  if no well-conditioned set can be found. Smoke test now 9/9.

### v2.5.17 addendum 2 — gate wired into the driver/job
* (A) SETUP-TIME pre-check: fep_driver validates the Boresch geometry
  (boresch._precheck_boresch_dict: r in [4,13] A, key angles in [30,150] deg)
  right after boresch.json is written and BEFORE the correction is finalized;
  a bad pick fails the run before any window is submitted.
* (B) RUNTIME gate: build_lsf_array now emits a post-equilibration geometry
  check into every complex-leg .lsf -- AFTER the eq temperature gate, BEFORE
  production. It re-validates the SAME six atoms against eq.rst via
  boresch_autoselect.validate_masks and exits 72 if the equilibrated geometry
  drifted out of bounds. It is FAIL-FAST only (no on-node re-selection): the
  analytic standard-state correction is fixed at setup time, so swapping atoms
  mid-job would silently bias dG. setup_leg drops boresch_eqcheck.json (atom
  masks) so the gate can run; if masks are unavailable the gate self-skips.
* Generated .lsf passes `bash -n`; smoke test 9/9.

### v2.5.17 addendum 3 — runtime gate now covers GUI (--auto-boresch) runs
* fep_driver bridges the legacy dict-style Boresch atoms (1-based serials
  aA/bA/cA + A/B/C from select_boresch_atoms, used by the GUI ABFE-Amber path)
  to Amber @serial masks (lig=A/B/C, rec=aA/bA/cA) and persists them in
  boresch.json. setup_leg then emits boresch_eqcheck.json, so the in-job
  post-equilibration gate (B) activates for GUI runs too -- not just
  mask-based auto-selection. Verified @serial masks reproduce the dict's own
  r/alpha/theta geometry. Purely additive; never overwrites existing masks.

### v2.5.17 addendum 4 — runtime gate import-robustness + honest SKIP
* FIX: the in-job gate ran a bare interpreter and failed `import amber_md`
  ('No module named amber_md'), then a trailing echo printed a FALSE
  'equilibrated geometry OK'. The gate now (a) injects sys.path from
  AMBER_MD_HOME, PYTHONPATH, and a walk-up from $WD that locates the
  amber_md/ package, so the import works on the compute node; and (b) uses
  distinct exit codes -- 0=PASS (prints OK), 72=FAIL (stops production),
  3=SKIP (prints 'SKIPPED (NOT validated)', never a false OK). Verified by
  executing the heredoc as an isolated interpreter from inside $WD for all
  four outcomes (PASS via walk-up, PASS via AMBER_MD_HOME, FAIL, SKIP).
