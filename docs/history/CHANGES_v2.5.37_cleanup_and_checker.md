# v2.5.37 - dual-copy cleanup + restraint-leg reliability checker

Date: 2026-06-20

(b) Removed amber_md/abfe_restraint_topology.py and tools/restraint_dualcopy_smoketest.py (dead since Option A in v2.5.36). The fep_driver dual-copy build block is gone; build_restraint_topology=True is now refused with an error rather than silently mis-coupling. SOURCES.txt + README_CLEAN_RUN updated.

(a) New tools/check_restraint_leg.py: parses dvdl_summary.csv + dG_estimators.csv, prints the dV/dl profile, and enforces the v2.5.36 gate (max|dV/dl|<=200, |BAR-TI|<=50). Exit 0=reliable, 1=unreliable, 2=IO error.
