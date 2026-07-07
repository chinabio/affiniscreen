CHANGES_rbfe_amber_planner (v2.5.3 / final57)
=============================================
Date: 2026-06-09

WHAT
----
Enable the Amber RBFE path in the GUI (Option A: planner stage), replacing the
hard-disabled "RBFE/Amber is a stub" gate.

BACKGROUND
----------
amber_md/rbfe_map.py is a complete, tested RBFE perturbation-map generator
(OpenFE + RDKit backends; redundant/mst/star/lomap networks; MCS atom mapping;
cost/score edge selection). It already writes nodes.csv, edges.csv,
mapping.json, map.graphml, diagnostics.txt and an executable run_edges.sh.
It was simply never wired into the Setup & Launch wizard, which advertised
RBFE/Amber as a stub and pointed users to OpenFE.

CHANGES
-------
amber_md/gui/pages/0_Setup_and_Launch.py
  1) COMPAT[("RBFE","Amber")] flipped False -> True with a descriptive note.
  2) New _build_commands() branch for (RBFE, Amber): runs
       python -m amber_md.rbfe_map -i <ligands> -o <wd>/rbfe_map
         --backend {auto|openfe|rdkit} --map-type {redundant|mst|star|lomap}
         --min-degree N --resname-a :L1 --resname-b :L2 [--hub H] [--use-hs]
     Params are read from the wizard (rbfe_backend, rbfe_map_type,
     rbfe_min_degree, rbfe_hub, rbfe_use_hs, rbfe_resname_a/b) with safe
     defaults, so it works out of the box even before dedicated widgets exist.

SCOPE / LIMITATION (read this)
------------------------------
This enables the PLANNER stage only. The generated run_edges.sh calls
  fep_driver --mode legacy ... --submit
per edge but ASSUMES the per-edge dual-topology prmtops (complex/solvent)
already exist. Automating that per-edge dual-topology build (tleap/pmx/femto)
is a separate, larger effort -- the future "--mode rbfe" driver (Option B).
Until then, RBFE/Amber yields a validated network + ready-to-edit run_edges.sh;
it does not by itself queue alchemical MD.

No change to MM-GBSA, ABFE, or any OpenFE path.

FILES CHANGED
-------------
- amber_md/gui/pages/0_Setup_and_Launch.py
