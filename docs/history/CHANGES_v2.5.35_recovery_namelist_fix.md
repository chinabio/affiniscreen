# v2.5.35 - recovery namelist fix + durable log

Date: 2026-06-20

## Bug (lig_12944901 complex_restraint lambda=0.800, job 4699289)
v2.5.34 built regrid/settle inputs by sed-editing prod.in -> 'tempi=, temp0='
(empty) -> pmemd "Cannot match namelist object name" -> instant crash. All 10
regrids AND the CPU settle fired (artifacts present) but every one died at parse.

## Fixes
1. _write_recovery_mdin(): parse clambda/temp0/timask/ifsc/nmropt from prod.in once,
   emit a clean self-contained NPT namelist. No more empty values / dangling commas.
2. mpirun --mca btl_tcp_if_include lo (kills em1.1720 warning that aborted MPI).
3. Durable recovery.log via _rec()/tee in every lambda dir (LSF .err was truncated).
