CHANGES_lock_openmm_mmgbsa_engine_marker (v2.5.7 / final63)
==========================================================
Date: 2026-06-10

CONTEXT
-------
final62 inferred MM-GBSA engine from GUESSED markers (mmgbsa/openmm.json,
openmm_system.xml) the OpenMM runner never writes. Audit of
amber_md/mmgbsa_openmm.py: OpenMM-only artifacts are work/production.nc and
work/md_log.txt; the per-run 'plan' dict was logged but NOT persisted. The
Amber pmemd path (batch.py/fep.py) writes neither.

FIX
---
amber_md/mmgbsa_openmm.py
  * import json (was missing).
  * After scoring, write self-describing mmgbsa/engine.json:
    {"engine":"OpenMM","method":"MM-GBSA","runner":"amber_md.mmgbsa_openmm",
     ...plan...}, next to FINAL_RESULTS_MMPBSA.dat. Best-effort.
amber_md/gui/results_lib.py
  * NEW mmgbsa_engine(wd) -> "OpenMM"|"Amber": (1) engine.json marker,
    (2) production.nc / md_log.txt (pre-final63), (3) default Amber.
amber_md/gui/pages/5_Results_Compare.py
  * Use rl.mmgbsa_engine(wd) instead of the guessed markers.

IMPACT
------
* New OpenMM MM-GBSA runs labelled via a real, self-describing marker;
  pre-final63 OpenMM runs still detected via production.nc/md_log.txt.
* Amber MM-GBSA -> "Amber". Detector centralized in results_lib for reuse.

FILES
-----
- amber_md/mmgbsa_openmm.py
- amber_md/gui/results_lib.py
- amber_md/gui/pages/5_Results_Compare.py
