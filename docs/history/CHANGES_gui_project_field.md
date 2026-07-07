# v2.5.0 build final30: editable LSF project field in the GUI

Trigger: "There is no way to change the project name on the GUI. Add a place to
default as your-project but user can change it."

## What changed
gui/pages/0_Setup_and_Launch.py
  * NEW widget in the "HPC / scheduler" expander:
        LSF project (-P)  [your-project]   (key=wiz_project)
    Editable; defaults to the lab allocation. Stored as params["project"].
  * GPU-queue submit command now appends `--project <value>` (stripped; omitted
    if the user blanks it).
  * OpenFE RBFE path already reads P.get("project", "your-project")
    (added final29), so the same field now drives it.

gui/openfe_common.py
  * Hardened the -P guard to strip whitespace, so a blanked field can never emit
    `#BSUB -P    ` (only a real, non-empty project produces a -P line).

## Flow (single field -> every path)
  HPC expander "LSF project (-P)" -> params["project"]
     -> amber MM-GBSA / ABFE GPU cmd:  --project <value>
     -> OpenFE RBFE:                   OpenFESettings(project=<value>) -> #BSUB -P
     -> mmgbsa_openmm / fep_driver:    consume --project (else HPCConfig default)

## Verified
  default            -> --project your-project ; openfe -P set
  user override      -> --project <theirs>              ; openfe -P set
  blanked            -> --project omitted               ; openfe -P omitted (strip)
  key missing (old)  -> falls back to your-project
  All .py compile.
