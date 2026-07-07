# v2.5.31d -- prod box-drift root cause: Boresch reference vs coords mismatch
eq CLEAN (T=298.5, rho=0.994, P~0). prod step-1 abort x5. boresch.RST r0=5.135A but
actual L1-P1 ~65A -> RESTRAINT 50,022 kcal/mol -> GPU halt. Fixes: (A) write-time
_verify_boresch_reference (min-image, fail loud), (B) prod vlimit, (C) GPU regrid
retry (no CPU), (D) gate hard-stop when amber_md missing on node. All on GPU.
