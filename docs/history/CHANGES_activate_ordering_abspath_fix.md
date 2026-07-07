# CHANGES - v2.5.3: activate-script ordering & absolute-path fix

## Symptom
A GPU MD job completed MD but the automatic MM/GBSA analysis never ran. The
job's .err contained:

    .../<job>.shell: line 17: ./activate_amber_md.sh: No such file or directory

MD still finished (pmemd.cuda comes from `module load amber/22.8`), but the
`amber-md` env required by MMPBSA.py was never activated, so no
FINAL_RESULTS_MMPBSA.dat was produced.

## Root cause (two compounding bugs)
1. `set -euo pipefail` was emitted AFTER `source <venv_activate>`, so the
   failed source was ignored and the script continued.
2. `config.HPCConfig.venv_activate` defaults to a RELATIVE
   "./activate_amber_md.sh"; LSF runs from ~/.lsbatch/... so it is not
   found. run_amber.py rewrote it, but GUI/batch paths did not.

## Fix
submit.py._header() now (a) emits strict mode first, (b) resolves
venv_activate to an absolute path, (c) guards its existence with a clear
error. config.py annotates the retained relative default.

## Impact
No need to re-run completed MD. New submissions activate correctly and fail
loudly if the activate script is ever missing.
