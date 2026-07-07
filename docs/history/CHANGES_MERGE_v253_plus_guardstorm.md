# MERGE: v2.5.3 (uploaded: submit.py activate/abspath fix) + final47 guard fix

## Base
User-uploaded amber_md_workflow_v2.5.3_patched.zip = my final46
(persisted-resname reload + fail-fast 0-atom mask guard in fep_driver) PLUS a
v2.5.3 submit.py patch:
  * _header() now emits `set -euo pipefail` BEFORE module load/source.
  * venv_activate resolved to an ABSOLUTE path with an existence guard
    (relative ./activate_*.sh was never found from ~/.lsbatch -> silent env
    failure that left MM/GBSA without its conda env). Touches config.py +
    submit.py only.

## Overlaid (was MISSING from the upload)
* tools/abfe_resume_guard.sh -> v3 (final47): SUBMISSION-AWARE. Skips a sweep
  while this run still has PEND/RUN jobs (read from fep/job_ids.json), so it can
  never again stack duplicate submissions (the 889-job storm).
* CHANGES_guard_submission_storm.md.

## Why the merge is clean
Disjoint files: the v2.5.3 patch is config.py + submit.py; the guard fix is a
shell script. fep_driver.py is byte-identical to final47 (all final46 driver
hardening intact). No conflicts.

## Net result = everything:
  final38-45 analysis/GUI + final46 resume-resname & 0-atom guard
  + v2.5.3 activate/abspath env fix + final47 guard anti-storm.
