# v2.5.27 -- HOTFIX: missing -ref on restrained stages (regression from 2.5.25)

SYMPTOM: complex_restraint failed immediately. heat stage (lambda=0.000) aborted:
  Unit 10 Error on OPEN: refc
  STOP PMEMD Terminated Abnormally!  (at "5. REFERENCE ATOM COORDINATES")

CAUSE: v2.5.25 turned on positional restraints (ntr=1) for heat and (by default)
dens+eq, but the run_stage commands were never given a -ref reference-coordinate
file. With ntr=1 pmemd MUST open a reference (default name 'refc'); none existed.
min was unaffected (ntr=0). This would have hit every leg's heat/dens/eq.

FIX: each restrained stage now passes -ref equal to its own starting coords:
  heat  -c min.rst   -ref min.rst
  dens  -c $DENS_C   -ref $DENS_C       ($DENS_C = heat.rst, or min.rst if no heat)
  eq    -c dens.rst  -ref dens.rst
Applied in BOTH the LSF-array and HREMD/resume run scripts (6 lines total).
Restraints hold atoms at the position they START each stage -- standard practice.

VALIDATION (no Amber): fep.py parses clean; exactly 6 '-ref' tokens; every
heat/dens/eq run_stage line carries -ref; min carries none.
RECOMMENDED: re-extract ONE window (lambda_0.000) and run min->heat->dens->eq->prod
before launching the full grid.
