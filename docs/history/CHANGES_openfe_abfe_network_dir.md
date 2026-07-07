CHANGES_openfe_abfe_network_dir (v2.5.4 / final60)
==================================================
Date: 2026-06-10

SYMPTOM
-------
Running an ABFE/OpenFE campaign submitted NOTHING (for any molecule). The FEP
Campaign page reported "No transformations under network_setup", even though
the ABFE planner had successfully written one transformation per ligand.

ROOT CAUSE
----------
ABFE/OpenFE plans into  <work_dir>/abfe_setup/transformations/*.json
RBFE/OpenFE plans into  <work_dir>/network_setup/transformations/*.json
(see 0_Setup_and_Launch.py: ABFE --out abfe_setup; RBFE --out network_setup).

But the campaign machinery hard-coded "network_setup":
  * openfe_campaign.CampaignLayout.at() -> network_dir = wd/"network_setup"
  * openfe_campaign.find_mapping_artifacts() search_dirs
  * 1_Job_Monitor.py status snapshot
So for ABFE, edges()/list_transformations() globbed an empty/nonexistent dir ->
task_list was empty -> submit_campaign() bsub'd nothing. The submitter itself
was correct (it iterates ALL transformations x replicates); it was simply
pointed at the wrong directory.

FIX
---
openfe_campaign.py
  * NEW resolve_network_dir(work_dir): returns whichever of network_setup/ or
    abfe_setup/ actually contains transformations (prefers a populated
    transformations/*.json; falls back to network_setup for back-compat).
  * CampaignLayout.at() now uses resolve_network_dir().
  * find_mapping_artifacts() search_dirs now include abfe_setup.
1_Job_Monitor.py
  * status snapshot uses resolve_network_dir() so ABFE edges are visible.

IMPACT
------
* ABFE/OpenFE campaigns now find all per-ligand transformations and submit one
  quickrun job per ligand (x n_replicates), as intended.
* RBFE/OpenFE unchanged (network_setup still resolved first).
* Pure GUI/orchestration change; no planner or protocol code touched.

FILES CHANGED
-------------
- amber_md/gui/openfe_campaign.py
- amber_md/gui/pages/1_Job_Monitor.py
