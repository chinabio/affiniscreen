# Release v2.5.61 — supporting documents

This folder bundles the investigation and patch artifacts behind the v2.5.61
timestep fix. The fix itself is ALREADY APPLIED in this package; the scripts here
are kept for provenance / manual re-application on other installs.

| File | What it is |
|---|---|
| FEP_dt002_crash_investigation_HANDOFF.md | Full crash-investigation handoff: dt=0.002 integrator instability on lig_12944901, multi-seed evidence, the GUI dt/nstlim decoupling, and the decision trail. |
| apply_dt001_GUI_abfe_fix.sh | Standalone patch for the GUI ABFE/Amber launch path (config.py dt_ps/prod_dt_ps + 0_Setup_and_Launch.py nstlim /2->/1 + fep_common widget). Verified against amber_md v2.5.60. |
| apply_dt001_clean_run_fix.sh | Earlier standalone patch targeting the run_amber.py / clean-run path. Superseded by the GUI script for GUI launches; kept for reference. |

NOTE: you do NOT need to run these scripts on a v2.5.61 install — the edits are
baked in. See ../../CHANGES_v2.5.61_timestep_and_restraint_analysis.md.
