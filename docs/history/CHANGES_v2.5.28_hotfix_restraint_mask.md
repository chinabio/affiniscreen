# v2.5.28 -- HOTFIX: restraint mask wildcard + duplicate heat nmropt/&wt

SYMPTOM (after v2.5.27 -ref fix): heat lambda=0.000 aborted:
  Error in group input::atommask.f::residue_namelist  /  unknown symbol:*
CAUSE: pmemd legacy group-input mask parser (ntr=1 restraintmask) rejects '*'.
FIX: posres_mask_default = "!:WAT,HOH,Na+,Cl-,K+,Mg2+,Ca2+,Zn2+ & !@H=" (no '*').
VALIDATED with parmed 4.3.1 vs complex_restraint.parm7 (272,652 atoms): broken &
fixed masks select the SAME 6,850 atoms; :1/:2=LIG/LIG (37 each, matched); Boresch
atoms on LIG (C5/F1/C16)+protein (ARG693/LEU695/MET699). Earlier 'complex.prmtop'
(148,238, res1=MET) was the plain topology, not the leg input -- symlink confirms
the leg uses complex_restraint.parm7. Dual-copy was never the bug.
ALSO rewrote _heat_in -> single nmropt=1 + single &wt END.
NEXT: re-extract lambda_0.000, run min->heat->dens->eq->prod, expect clean exit.
