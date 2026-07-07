# v2.5.29 -- HOTFIX: malformed heat &wt TEMP_0 ramp broke heat on ALL legs

SYMPTOM: all five legs failed at heat, lambda=0.000, before any dynamics:
  WEIGHT CHANGES:
  Error: Invalid TYPE flag in line:
  TEMP_0        0  80000    5.000000  298.000000      0      0
  STOP PMEMD Terminated Abnormally!
REGRESSION: decharge/vdw ran fine before the v2.5.25 heat stage existed.

ROOT CAUSE: pmemd's nmropt=1 weight-change reader needs the &wt namelist OPENER on
its own line. The heat template emitted it inline:
  &wt type='TEMP_0', istep1=0, ... /     -> mis-tokenized -> "Invalid TYPE flag".

FIX (heat stage only): use the validated FEP-SPell-ABFE layout already present in
amber_md/equilibration_fepspell.py (lines 86-93):
   /
   &wt
   type = 'TEMP_0',
   istep1 = 0, istep2 = <ramp>,
   value1 = <Ti>, value2 = <T0>,
   /
   &wt
   type = 'END',
   /
  DISANG=boresch.RST    (only when Boresch active)
The dens/eq/prod END-only &wt tail is UNCHANGED (proven-working for years).

VALIDATED: heat.in rendered in-sandbox; &wt opener on its own line; fep.py parses.
The fix chain to date: 2.5.26 version track, 2.5.27 -ref, 2.5.28 mask wildcard,
2.5.29 &wt card -- each uncovered the next latent defect in the new heat stage.
This is the last heat setup-time bug; run ONE window before the full grid.
