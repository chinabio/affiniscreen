# v2.5.40 - heat-stage Option A hotfix

Date: 2026-06-20

Real-run failure (lig_12944901 complex_restraint, lambda=0.000): min PASSED but heat aborted:
  TI Mask 1 :1; matches 19 atoms / TI Mask 2 :2; matches 24 atoms
  ERROR: timask1/2 must match the same number of atoms for non-softcore run

Root cause: _heat_in() has its OWN mdin template; v2.5.36 routed min/dens/eq/prod through _ti_kw_line but missed heat, so heat still wrote icfe=1 + timask=':1'/':2'. Those masks are meaningless on the real single-copy complex.prmtop (residues 1/2 are protein), hence the abort.

Fix: _heat_in now uses _ti_kw_line(cl), which returns '' for the restraint stage. Confirmed: system.prmtop correctly points at build/complex.prmtop (single copy) and boresch.RST already carries real-topology indices with rk2=rk3=0 at lambda=0 -- Option A topology routing was correct; only the heat mdin needed the fix.
