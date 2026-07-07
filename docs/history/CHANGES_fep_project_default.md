# v2.5.0 build final28: fix FEP default LSF project

## Bug
amber_md/fep_driver.py defined:
    p.add_argument("--project", default="default")
and passed it straight into HPCConfig(project=a.project, ...). So an Amber FEP
submission WITHOUT an explicit --project emitted:
    #BSUB -P default
i.e. the bogus project literal "default" instead of the real allocation
your-project. (The MM-GBSA path was already correct: --project default
None, falling back to HPCConfig.project.)

## Fix
* CLI default changed to None.
* HPCConfig is now built WITHOUT project unless --project was explicitly given,
  so the dataclass default (your-project) wins by default; an explicit
  --project still overrides it.

## Verified
  _build_configs(project=None)      -> hpc.project == your-project
  _build_configs(project=explicit)  -> hpc.project == explicit
  (never "default")

## Reminder: all interactive/bsub commands need the project too
    bsub -P your-project -q gpu -gpu "num=1" -W 0:15 -Is bash
