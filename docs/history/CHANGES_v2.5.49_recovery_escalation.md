# v2.5.49 - escalating box-drift recovery

Date: 2026-06-21

Run abfe_20260620_222307 lig_12944901 complex_restraint: 15/16 windows COMPLETED (recovery works: lambda 0.300 took 5 attempts, 0.900 took 3). Only lambda=0.150 FAILED -- a persistent steric clash (VDWAALS +20,000, PRESS 8762, TEMP 361 at step 19000) that 10 identical 50 ps @ 1 fs regrids never resolved because each regrid re-seeded near the same strained geometry.

Changes:
  1. Escalating regrid ladder keyed on a drift counter:
       tier1 (drift 1-2): 50 ps @ 1.0 fs, taup=5, from good_rst
       tier2 (drift 3-4): 100 ps @ 0.5 fs, taup=2, from orig_eq.rst (pristine eq)
       tier3 (drift >=5): CPU density-settle, then resume GPU
  2. orig_eq.rst preserves the TRUE equilibrated restart (the prior code overwrote the file named eq.rst on every regrid, permanently losing real eq coords).
  3. Clearer recovery.log: tier + step/dt + actual source file (old log always said 'from eq.rst' even though good_rst had been overwritten).

Caveat: a real clash may still need finer lambda spacing (insert 0.10/0.20). Re-run ONLY lambda=0.150 -- the other 15 prod.out are complete.
