# v2.5.25 -- decouple/vdw reliability: NVT heat + restrained eq + vlimit
New stage order min->heat->dens->eq->prod. Soft-core block already matched
FEP-SPell-ABFE exactly; gaps fixed were staging (no NVT heat), unrestrained eq,
and no vlimit valve. See README_CLEAN_RUN.md / v2.5.26 changes.
