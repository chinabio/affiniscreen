# v2.5.43 - driver-level restraint-input guard (RBFE + ABFE)

Date: 2026-06-20

setup_leg() already raises on a TI-keyword leak (v2.5.42). v2.5.43 adds an explicit scan in run_fep right after setup_leg for stage=='restraint', logging 'restraint-input guard: OK' or aborting with return 2 and a per-leg list of offenders. RBFE coverage is automatic since each edge runs through run_fep via _run_single_edge.
