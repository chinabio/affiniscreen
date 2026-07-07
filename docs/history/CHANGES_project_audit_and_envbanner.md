# v2.5.0 build final29: cluster-wide project audit + preflight env banner

Trigger: "any bsub needs project your-project."

## Submit-path project audit (every bsub/#BSUB -P emitter)
  fep.py            -> #BSUB -P {hpc_cfg.project}        OK (uses HPCConfig)
  fep_driver.py     -> None -> HPCConfig.project         OK (fixed final28)
  mmgbsa_openmm.py  -> a.project or hpc.project; guarded OK (+warn added)
  submit.py         -> #BSUB -P {cfg.project}            OK
  run_amber.py      -> None -> HPCConfig conditional      OK
  GUI 2_MMGBSA / 3_FEP / openfe_common text_input default "your-project" OK

### Bugs found & fixed
1. gui/pages/0_Setup_and_Launch.py (OpenFE RBFE path) hardcoded project="" in
   OpenFESettings, so the generated OpenFE LSF script emitted a malformed
   `#BSUB -P ` (empty). FIX: project=P.get("project","your-project").
2. gui/openfe_common.py emitted `#BSUB -P {s.project}` UNCONDITIONALLY. FIX:
   only emit the -P line when project is non-empty; warn (and skip) otherwise,
   so an empty project can never produce a broken header.
3. mmgbsa_openmm.py: added a warning when no project resolves (it already
   correctly omitted the -P line via `if project:`).

## Preflight env banner (preflight_abfe.py)
Prints + records in abfe_preflight_report.json:
  CONDA_DEFAULT_ENV, AMBERHOME, which pmemd.cuda, LOADEDMODULES
and notes if the conda env is not 'amber-md'. So the report self-documents that
it ran in the right place (amber-md + amber/22.8).

## Verified
* openfe header: project set -> `#BSUB -P your-project`; empty -> no -P
  line at all (no malformed blank). 
* preflight env banner prints and lands in JSON.
* All .py compile.

## Reminder
Interactive shells still need it explicitly:
    bsub -P your-project -q gpu -gpu "num=1" -W 0:15 -Is bash
