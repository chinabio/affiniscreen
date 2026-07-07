# v2.5.0 build final37: MM-GBSA failure diagnostics

## Why
A user MM-GBSA run (mmgbsa_20260607_131542, predates the resname change) failed
with a GENERIC CalcError: 'mmpbsa_py_energy failed with prmtop complex.prmtop'
on the COMPLEX GB contribution, across all MPI ranks (MPI_ABORT). The real
cause lives in the per-calc temp file (_MMPBSA_complex_gb.mdout.0), but
keep_files=0 + MPI cleanup removed them -> 'ls _MMPBSA_*' found nothing, so the
root cause could not be read.

## Fix (amber_md/mmgbsa.py)
* write_mmgbsa_input: keep_files=0 -> keep_files=1 so a FAILED run retains the
  per-calculation mdout/out files for inspection.
* Self-diagnose block hardened: on non-zero MMPBSA.py exit it now
    - globs _MMPBSA_*.out*, _MMPBSA_*.mdout*, reference.frc,
    - COPIES them into <workdir>/mmpbsa_failed_tmp/ BEFORE MPI cleanup can
      delete them,
    - logs that directory + tails the first 6 files into the .log.

This does NOT change the science; it only preserves and surfaces the real error.

## Next run will show, in the .log and in mmpbsa_failed_tmp/:
  the actual sander/mmpbsa_py_energy error for the complex (radii / charge /
  atom-count / mask mismatch), which tells us the true root cause.
