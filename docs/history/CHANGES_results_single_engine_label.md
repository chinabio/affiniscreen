CHANGES_results_single_engine_label (v2.5.8 / final64)
=====================================================
Date: 2026-06-10

amber_md/gui/pages/4_Results_Single.py
  * MM-GBSA headline now shows the detected ENGINE alongside DG_bind, using the
    locked detector rl.mmgbsa_engine(wd) (final63): engine.json marker ->
    OpenMM MD artifacts (production.nc/md_log.txt) -> default Amber.
  * Adds an "Engine" metric + a Method/Engine caption. No other behaviour
    changed; OpenFE and Amber-FEP headlines untouched.

FILES
-----
- amber_md/gui/pages/4_Results_Single.py
