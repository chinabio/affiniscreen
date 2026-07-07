CHANGES_rbfe_amber_endtoend (v2.5.3 / final58)
==============================================
Date: 2026-06-09

Option B: full end-to-end Amber RBFE. Closes the final57 gap (planner only):
per-edge dual-topology prmtops are now built automatically and the TI legs are
submitted, so RBFE/Amber is queue-able like ABFE.

NEW  amber_md/rbfe_topology.py
  build_rbfe_edge_topology(): for edge A->B, parametrises BOTH ligands as :L1
  and :L2 and builds two solvated dual-topology systems (complex: protein+L1+L2;
  solvent: L1+L2). Reuses PDBCleaner/LigandParametrizer/protonation/_WATER_MAP/
  SystemConfig. Masks always :L1/:L2 (match rbfe_map).

CHG  amber_md/fep_driver.py
  --mode gains "rbfe"; new args --ligands/--rbfe-map-dir/--edges-csv/--only-edges.
  run_rbfe() reads edges.csv and per edge: extract A,B records -> build dual
  topology -> delegate to run_fep() LEGACY (runs the dual-topology TI edge and
  computes relative_binding_dG = dG_complex - dG_solvent with --analyze).
  Writes RBFE_CAMPAIGN.json.

CHG  amber_md/gui/pages/0_Setup_and_Launch.py
  COMPAT note -> "end-to-end". (RBFE,Amber) branch now emits rbfe_map THEN
  fep_driver --mode rbfe (with --analyze; --submit on GPU queue). Launch UI
  recognises RBFE/Amber+GPU queue as an LSF submission.

THERMO  per edge ddG(A->B)=dG_complex(L1->L2)-dG_solvent(L1->L2); network cycle
  closure handled downstream by the results pages/aggregator.

LIMITATIONS  ligands placed at own coords (not yet MCS-RMSD-fit pre-tleap;
  mapping.json gives correspondence for a future refinement); minimisation
  relaxes initial overlap. No change to MM-GBSA/ABFE/OpenFE paths.
