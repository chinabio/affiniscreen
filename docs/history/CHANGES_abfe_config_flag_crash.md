CHANGES_abfe_config_flag_crash (v2.5.3 / final55)
=================================================
Date: 2026-06-09

SYMPTOM
-------
Launching an Amber ABFE job from the GUI (Setup & Launch) died immediately
with, in launch_abfe_*.log:

    amber-fep: error: unrecognized arguments: --config .../wizard_config.json

No alchemical window was ever submitted.

ROOT CAUSE
----------
amber_md/gui/pages/0_Setup_and_Launch.py::_build_commands() built the Amber
ABFE command with `--config <wizard_config.json>`. The fep_driver entrypoint
(amber_md/fep_driver.py) has NO --config option; it derives the ABFE alchemical
masks from --ligand-resname and reads all other settings from explicit CLI
flags. argparse therefore rejected the unknown flag and exited before run_fep().

A latent second crash was lurking in the same block: it appended
_protonation_flags(P) (--no-protonation / --protonate ...), which fep_driver
also does not accept. Any user who toggled protonation off, or added an
override, would have hit the same "unrecognized arguments" abort.

A third, cosmetic-but-misleading issue: the written wizard_config.json carried
FEPConfig's RBFE-style default masks (timask1=:LIA, timask2=:LIB). Those
residues do not exist in an absolute (single-topology) ABFE prmtop. The driver
ignored them (it re-derives masks), but anyone reading the config would be
misled.

FIX
---
1) 0_Setup_and_Launch.py (ABFE/Amber branch):
   - Removed `--config <cfg_path>` from the fep_driver command.
   - Removed `c += _protonation_flags(P)` from the ABFE command.
   - Still write wizard_config.json as a PROVENANCE record (not passed to the
     driver), now with the correct ABFE ligand resname.
   - Pass the system-prep knobs the driver DOES accept so GUI choices are not
     silently dropped: --ligand-charge, --charge-method, --box-buffer,
     --project.
   (The MM-GBSA / run_amber.py branch still uses --config; that entrypoint
    legitimately supports it. Unchanged.)

2) amber_config.py (build_config):
   - For method == "ABFE", set single-topology masks in the provenance config
     to match fep_driver._derive_abfe_masks(): timask1/scmask1 = :<resname>,
     timask2/scmask2 = "", crgmask = :<resname>. RBFE path unchanged.

IMPACT
------
- Amber ABFE now launches from the GUI without the --config crash.
- No behavioural change to MM-GBSA, RBFE, or the OpenFE paths.
- The masks/charge/buffer/project the GUI shows now actually reach the run.

FILES CHANGED
-------------
- amber_md/gui/pages/0_Setup_and_Launch.py
- amber_md/gui/amber_config.py
