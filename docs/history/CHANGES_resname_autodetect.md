# v2.5.0 build final34: ligand resname auto-detect + ABFE mask derivation

## Problem
Resname was effectively hardcoded. The smoke-test set --ligand-resname but the
alchemical masks --timask1/--scmask1 stayed at :LIG, so a 'UNK' ligand (the
user's actual mol2) would decouple a NON-EXISTENT residue -> garbage/failed ABFE.

## Fix (fep_driver.py)
* detect_ligand_resname(ligand_file): reads the residue-name column of a .mol2
  @<TRIPOS>ATOM block (e.g. 'UNK'); .sdf/.mol -> default. Never raises.
* --ligand-resname now defaults to None. Resolution priority in run_fep (abfe):
    explicit --ligand-resname  >  auto-detect from .mol2  >  'LIG'
* _derive_abfe_masks(a): for ABFE, derives ALL masks from the resolved resname
  UNLESS the user explicitly set one on the CLI:
    timask1=scmask1=crgmask=:<resname>,  timask2=scmask2=''   (single-topology)
  crgmask is now threaded into FEPConfig.
* Logs the resolved resname (with source) and the final masks.

## Fix (run_abfe_smoketest.sh)
* RESNAME is now OPTIONAL (3rd arg). If omitted, the driver auto-detects.
* Removed the hardcoded --ligand-resname/mask coupling; passes
  --ligand-resname only when the user supplies one.

## Verified (against the user's real UNK mol2)
  detect mol2 -> ('UNK','mol2');  sdf -> ('LIG','default')
  A auto      -> resname UNK, timask1/scmask1/crgmask=:UNK, timask2/scmask2=''
  B override  -> :MYL everywhere
  C sdf       -> LIG default
  D user mask -> explicit --timask1 respected; others still derived
  All .py compile.

## So for the user's ligand:
  bash run_abfe_smoketest.sh ~/Run_dir/protein.pdb ~/Run_dir/ligand.mol2
  # resname UNK auto-detected; masks become :UNK automatically.
