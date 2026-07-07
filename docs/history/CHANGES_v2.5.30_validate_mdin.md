# v2.5.30 -- tools/validate_mdin.py : catch mdin parser bugs BEFORE the GPU

WHY: v2.5.27, 2.5.28 and 2.5.29 each fixed a real defect in the new heat stage,
but each only surfaced the NEXT one, because pmemd aborts on the FIRST setup error
it meets -- so every fix cost another cluster round-trip. This validator parses the
generated .in files locally and reports ALL such issues at once, in milliseconds.

CHECKS:
  1. &wt cards          -- the '&wt' opener must be alone on its line (else pmemd
                           'Invalid TYPE flag'); nmropt=1 needs an END terminator.   [2.5.29]
  2. restraintmask      -- no '*' residue wildcards (pmemd group parser rejects).      [2.5.28]
  3. -ref wiring        -- any ntr=1 stage's pmemd command must pass -ref, else
                           'Unit 10 Error on OPEN: refc'.                              [2.5.27]
  4. heat is NVT        -- ntp=0/ntb=1 (never barostat a heating soft-core ligand).
  5. icfe=1 TI masks    -- timask1/timask2 present and non-empty.
  6. DISANG             -- referenced restraint file exists in the window.

USAGE:
  python tools/validate_mdin.py /path/to/<leg>            # one leg
  python tools/validate_mdin.py /path/to/abfe_<timestamp> # whole run (all legs)
Exit code 0 = clean (safe to submit); non-zero = at least one fatal issue.

SELF-TEST (recorded): on a known-bad window it reports 4 errors (2 inline-&wt,
1 wildcard mask, 1 missing -ref); on a clean v2.5.29 window it reports 0 and PASS.

RECOMMENDED WORKFLOW: after the GUI generates a run, before launching:
  python tools/validate_mdin.py ~/Run_dir/run_v250/abfe_<ts>/lig_*/fep/*
then run ONE window (min->heat->dens->eq->prod) as a smoke test, then the full grid.
