# Architecture

The GUI is a thin orchestration layer over unchanged engine code. Selections on
**Setup & Launch** are translated into shell commands by `_build_commands()` and
spawned detached; nothing in the science layer is modified by the GUI.

## Method × Engine routing (GUI-supported)

```
MM-GBSA / Amber            -> run_amber.py (per-ligand fan-out)
MM-GBSA / OpenMM/OpenFE    -> amber_md.mmgbsa_openmm (per-ligand fan-out)
ABFE    / OpenMM/OpenFE    -> amber_md.abfe_openfe_plan  (OpenFE env python)
RBFE    / OpenMM/OpenFE    -> openfe plan-rbfe-network -> FEP Campaign run/gather/solve
```

Amber ABFE and Amber RBFE are disabled in the GUI's compatibility matrix
(`amber_md/gui/pages/0_Setup_and_Launch.py::COMPAT`). Their command builders now
raise a clear error if ever reached, and the FEP Campaign page no longer offers
an "Amber TI" engine. The engine modules themselves (`fep_driver`, `rbfe_map`,
`abfe_*`) remain importable and CLI-usable.

## Key modules

| Module | Role |
|--------|------|
| `amber_md/gui/Home.py` | Landing page + version badge. |
| `amber_md/gui/pages/0_Setup_and_Launch.py` | Unified wizard; `COMPAT` matrix; command builder. |
| `amber_md/gui/pages/3_FEP_Campaign.py` | OpenFE run/gather/solve/cycle-closure. |
| `amber_md/gui/openfe_campaign.py` | OpenFE campaign adapter. |
| `amber_md/gui/amber_campaign.py` | Amber campaign adapter (read-only result parsing; retained). |
| `amber_md/gui/routing_help.py` | Method×Engine×count routing table. |
| `amber_md/fep_driver.py` | FEP engine driver (CLI-supported ABFE/RBFE). |
