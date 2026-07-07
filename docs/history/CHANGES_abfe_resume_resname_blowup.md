# v2.5.0 build final46: fix ABFE whole-system-softcore blow-up on --resume

## Symptom
All 76 ABFE windows failed. dens.out: 45k+ lines "Molecule N is completely
softcore"; pmemd: "TI Mask 1 :LIG; matches 0 atoms" / "TI region 1: 148238".

## Root cause (RESUME bug, triggered by the guard wrapper)
mol2 resname is UNK. First launch -> masks :UNK, matched UNK prmtop (OK). The
guard's first --resume sweep ran fep_driver WITHOUT --ligand-file/-resname;
run_fep fell through to ligand_resname="LIG" and rewrote per-window masks to
:LIG. prmtop residue is UNK -> :LIG matches 0 atoms -> whole system softcore ->
all windows died.
  10:27 ligand resname: UNK ; masks :UNK  (OK)
  10:29 ligand resname: LIG (default; no ligand file); masks :LIG  (BROKE IT)

## Fixes
1. _load_persisted_resname(): reload ligand_resname from the prior build's
   abfe_topology_inputs.json before defaulting. --resume now reuses the exact
   prmtop residue name. Priority: explicit flag > persisted > mol2 > "LIG".
2. _assert_mask_matches()/_count_residue_atoms(): fail-fast guard after topology
   build; abort if ':<resname>' matches 0 atoms in complex/solvent prmtop
   instead of submitting 76 doomed windows. parmed or raw-prmtop fallback
   (unit-tested: UNK->11, LIG->0->ABORT).
3. tools/abfe_resume_guard.sh v2: resume sweep always passes --ligand-file +
   --ligand-resname; adds FIRST_SWEEP_DELAY_MIN (30) so first resume can't
   stampede a fresh launch.

## Note
abfe_production_20260608_102750 is unrecoverable (dens.in rewritten with :LIG).
Relaunch fresh; passing --ligand-resname LIG recommended.
