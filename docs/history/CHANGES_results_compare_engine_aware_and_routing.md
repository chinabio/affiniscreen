CHANGES_results_compare_engine_aware_and_routing (v2.5.6 / final62)
==================================================================
Date: 2026-06-10

(1) Results-Compare: engine/method awareness + RBFE-network rollup
  amber_md/gui/pages/5_Results_Compare.py
  * _collect() now emits an explicit ENGINE column and refined METHOD per row:
      - *_result.json -> method ABFE|RBFE, engine OpenFE
      - fep/ -> method ABFE|RBFE|FEP, engine Amber (+ 'detail' provenance)
      - mmgbsa/ -> method MM-GBSA, engine OpenMM (marker) else Amber
  * NEW RBFE-network rollup: detects a perturbation network and SOLVES per-edge
    ddG into per-ligand dG + cycle-closure QC inline, for BOTH engines:
      - Amber : <parent>/edges/<A~B>/fep/{complex,solvent}/summary.json
      - OpenFE: solved ddG TSV (openfe_campaign dg_tsv) if gathered
    Uses fep_common.solve_network / cycle_closure_residuals (matches Campaign).

(2) Workflow routing matrix help tab
  NEW amber_md/gui/routing_help.py
  * ROUTING_ROWS / routing_dataframe() (10 rows) + render_routing_help()
    expander, rendered near the top of Results-Compare; reusable anywhere.
  * Encodes the clarified routing incl. asymmetries (OpenMM MM-GBSA submits
    only from Setup & Launch; Amber ABFE multi fans out from Setup & Launch;
    RBFE network ranking lives on FEP Campaign).

No engine/driver/protocol code changed; GUI + one helper module.

FILES
-----
- NEW amber_md/gui/routing_help.py
- amber_md/gui/pages/5_Results_Compare.py
