# v2.5.36 Option A restraint leg + reliability gate

Restraint leg = single-copy plain MD (icfe=0) with lambda-scaled Boresch &rst (k=lambda*k_full). Replaces the v2.5.23 dual-copy TI whose :2 was an identical (non-decoupled) ligand copy -> 1/r^12 singularity (dV/dl -1209, BAR-TI ~500). build_restraint_topology defaults False. Analyzer gate: dG_reliable / reliability_reasons / max_abs_dvdl_kcal.
