CHANGES_multiligand_fanout (v2.5.4 / final59)
=============================================
Date: 2026-06-10

SYMPTOM
-------
With multiple molecules selected, Amber ABFE submitted only ONE job, for the
first molecule. (Amber MM-GBSA had the same latent flaw: one job for the whole
file, first record treated as the ligand.)

EXPECTED BEHAVIOUR
------------------
* MM-GBSA (any engine): one job PER molecule.
* ABFE   (any engine): one job PER molecule.
* RBFE   (any engine): ONE calculation using ALL molecules (a perturbation
  network / campaign). This was already correct and is unchanged.

ROOT CAUSE
----------
_build_commands() in 0_Setup_and_Launch.py:
  * ABFE/Amber passed the entire multi-record --ligand-file to a single
    fep_driver invocation with no --ligand-index; the driver defaults to
    index 0 and builds topology for the first record only. No loop -> 1 job.
  * MM-GBSA/Amber likewise passed the whole file to one run_amber.py call.
The OpenMM MM-GBSA branch already fanned out via split_ligand_file(); the Amber
branches simply never adopted that pattern.

FIX
---
0_Setup_and_Launch.py
  * ABFE/Amber: split the ligand file (split_ligand_file) and emit ONE
    fep_driver --mode abfe command per ligand into <wd>/lig_<name>/, each with
    its own single-record SDF and per-ligand provenance config.
  * MM-GBSA/Amber: same fan-out into <wd>/lig_<name>/, one run_amber.py per
    ligand.
  * _cmd_label(): also recognise the fep_driver spelling "--work-dir" (not just
    "--workdir") so per-ligand ABFE jobs get unique launch_<name>.log files
    instead of all writing launch_abfe_0.log.
  * RBFE/Amber + RBFE/OpenFE: UNCHANGED (correctly aggregate all molecules).

IMPACT
------
* Selecting N molecules now yields N MM-GBSA jobs or N ABFE jobs, each in its
  own lig_<name>/ directory (compatible with batch_aggregate ranking).
* RBFE still builds one network over all molecules.
* No driver/engine code changed; this is purely the GUI command builder.

FILES CHANGED
-------------
- amber_md/gui/pages/0_Setup_and_Launch.py
