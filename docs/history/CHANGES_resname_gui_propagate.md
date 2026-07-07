# v2.5.0 build final35: ligand resname auto-detect propagated to the GUI

Builds on final34 (driver-side auto-detect + ABFE mask derivation).

## Refactor: single shared detector
* Moved detect_ligand_resname() into core amber_md/utils.py (was a local copy
  in fep_driver). fep_driver now imports it from .utils -- ONE implementation.
* GUI Setup page imports the SAME amber_md.utils.detect_ligand_resname.
  Verified: fep_driver.detect_ligand_resname IS utils.detect_ligand_resname.

## GUI Setup page (0_Setup_and_Launch.py)
* New field "Ligand residue name (blank = auto)" (key wiz_resname), right under
  the ligand picker. Stored as params["ligand_resname"].
* LIVE preview caption:
    - .mol2 -> "Detected residue name in mol2: UNK -> effective: <eff>"
      (effective = user override if typed, else detected)
    - .sdf  -> "No residue name in this file (SDF) -- tleap will assign 'LIG'..."
* ABFE launch: passes --ligand-resname only when the field is non-blank;
  blank => the driver auto-detects (mol2) / defaults to LIG (sdf).

## How each method benefits
* ABFE/Amber : field -> --ligand-resname -> drives timask/scmask/crgmask
               (final34). Blank -> driver auto-detects from mol2.
* RBFE/OpenFE & MM-GBSA : these CLIs do not take a resname; the GUI preview
  shows the user WHICH residue tleap/antechamber will assign, so there are no
  silent surprises. (No behavioural change to those engines.)

## Verified
  utils detect: mol2->('UNK','mol2'), sdf->('LIG','default')
  driver shares the same function object
  GUI blank field -> no --ligand-resname (driver auto); set -> passes it
  GUI 'effective' caption: blank->UNK, 'MYL'->MYL
  All .py compile.
