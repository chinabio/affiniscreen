# v2.5.31 -- mdin validation wired into the driver (auto-gate before bsub)

amber_md/mdin_validator.py (importable core) + a GATE in fep_driver.run_fep that
validates each leg's freshly written .in files and run script BEFORE any bsub. On a
fatal issue the driver logs the exact cause and returns rc=2 without submitting.
tools/validate_mdin.py is now a thin CLI over the module. --skip-mdin-check overrides.
GUI Setup&Launch notes that validation runs automatically for Amber ABFE/RBFE.

Checks (each tied to a real past failure): inline &wt card (2.5.29), '*' restraintmask
wildcard (2.5.28), ntr=1 stage missing -ref (2.5.27), non-NVT heat, missing TI masks,
missing DISANG file.
