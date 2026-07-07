# v2.5.0 build final36: MM-GBSA enforces the ligand resname (auto-detected)

Builds on final35 (shared detector in utils; ABFE + GUI preview).

## What MM-GBSA already did
mmgbsa_openmm._force_ligand_resname() + antechamber -rn ALREADY force the
ligand residue name to MMGBSAConfig.ligand_resname. BUT run_amber.py
--lig-resname defaulted to "LIG" and was never auto-detected, so a 'UNK' mol2
was silently relabelled 'LIG' and the GUI never passed one.

## Fix
* run_amber.py: --lig-resname default "LIG" -> None; before WorkflowConfig,
  resolve  explicit > auto-detect(.mol2) > "LIG"  via the shared
  amber_md.utils.detect_ligand_resname. Resolved name is then ENFORCED
  downstream via antechamber -rn / PDB rewrite (unchanged).
* GUI Setup: MM-GBSA launch appends --lig-resname when the resname field
  (params['ligand_resname'], final35) is set; blank -> auto-detect.

## Net effect (all share ONE detector in amber_md.utils)
  ABFE/Amber  -> timask/scmask/crgmask (final34)
  MM-GBSA     -> antechamber -rn / PDB rewrite (THIS build)
  RBFE/OpenFE -> GUI preview only
  Blank everywhere -> auto-detect from .mol2 (else LIG).

## Verified (user's real UNK mol2)
  A mol2->UNK ; B override->MYL ; C sdf->LIG ; D no file->LIG
  GUI blank-> no flag ; GUI 'UNK'-> --lig-resname UNK ; all .py compile.
