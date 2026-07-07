CHANGES_fep_campaign_engine_selector (v2.5.5 / final61)
=======================================================
Date: 2026-06-10

The FEP Campaign page now drives BOTH engines behind one engine selector, so
OpenFE and Amber-TI share a single run -> gather -> solve -> cycle-closure UI.

NEW  amber_md/gui/amber_campaign.py
  Streamlit-free Amber-TI RBFE campaign adapter mirroring openfe_campaign's API
  (CampaignLayout.at, edges, task_list, campaign_status, submit_campaign,
  gather_command/gather_now, find_mapping_artifacts, read_edges_table,
  mapping_network_mermaid, preview_script, JOB_PREFIX, count_my_jobs).
  * Layout: <wd>/rbfe_map (planner) + <wd>/edges/<A~B> (run dirs).
  * edges(): from rbfe_map/edges.csv (falls back to existing edges/ dirs).
  * ddG(A->B) = dG_complex - dG_solvent, read from each edge's
    fep/complex/summary.json and fep/solvent/summary.json.
  * submit_campaign(): spawns ONE detached fep_driver --mode rbfe --analyze
    --submit (driver builds dual-topology + bsubs all per-edge legs). Needs
    protein+ligands (the files the map was planned from), from the sidebar.
  * gather_now(): in-process; writes <wd>/amber_rbfe_ddg.tsv.

CHG  amber_md/gui/pages/3_FEP_Campaign.py
  * Sidebar engine selector "OpenFE / OpenMM" | "Amber TI"; defaults to the
    last launch's engine. Dispatch via oc = oc_amber if Amber else oc_openfe.
  * Run tab branches (OpenFE throttle+per-edge bsub vs Amber single submit +
    protein/ligands inputs). Gather tab branches (openfe gather vs in-process).
  * Solve + cycle-closure tabs UNCHANGED (engine-agnostic fep_common consuming
    (a,b,ddG,err)). Status/Atom-mapping driven through the shared adapter API.

IMPACT
------
* Amber RBFE campaigns now run/monitor/solve from the same page as OpenFE,
  incl. per-ligand DeltaG ranking + cycle-closure QC.
* OpenFE path behaviourally unchanged. No engine/driver/protocol code changed.

NOTES
-----
* Amber per-edge ddG uncertainty is a 0.5 kcal/mol placeholder until per-leg
  stderr is propagated from FEPAnalyzer into each leg summary.json.
* Amber ABFE (per-ligand jobs) stays on Setup & Launch + Job Monitor; this page
  covers the network (RBFE) flow for Amber.

FILES
-----
- NEW amber_md/gui/amber_campaign.py
- amber_md/gui/pages/3_FEP_Campaign.py
