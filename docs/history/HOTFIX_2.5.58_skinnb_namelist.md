# HOTFIX 2.5.58 - skinnb must live in &ewald, not &cntrl

## Regression
v2.5.54 introduced `skinnb=3.0` to cure the spurious pmemd.cuda
"Periodic box dimensions have changed too much" halt. It was placed in the
**&cntrl** namelist. pmemd's &cntrl parser does not define `skinnb`, so AMBER 22
aborted immediately:

```
At line 975 of file .../pmemd/src/mdin_ctrl_dat.F90 (unit = 5, file = 'dens.in')
Fortran runtime error: Cannot match namelist object name skinnb
STAGE FAILED: leg=complex_restraint lambda=0.000 stage=dens rc=2
```

Every window in 2.5.54-2.5.57 died at the **dens** stage in ~2 min ->
`summary.json: n_windows=0, complete=false`. No MD ran.

## Root cause
`skinnb` is a **Particle-Mesh-Ewald variable**, defined in the **&ewald**
namelist (pmemd default `2.0`). Confirmed by:
- AMBER developers mailing list (the canonical report of this exact GPU error):
  the documented fix is *"increase skinnb to 3.0 ... by adding &ewald namelist"*.
- ParmEd's namelist dictionaries: `skinnb` is listed under `ewald` (pmemd),
  not `cntrl`.

## Fix
Emit a dedicated block after &cntrl in eq / dens / prod / restraint-prod mdins
and in the recovery mdin:

```
&cntrl
  ... (no skinnb here) ...
/
&ewald
  skinnb=3.0,
 /
```

`_min_in` is intentionally left unchanged (minimization has no PME pairlist
dynamics). All four production mdin renders were validated programmatically:
&cntrl and &ewald are separate, both namelists are closed, and `skinnb` appears
only in &ewald.

## Action required
**Upgrade to 2.5.58 before any run.** 2.5.54-2.5.57 cannot complete a single
MD step. The physics intent (larger pairlist skin to avoid the benign-NPT GPU
cell-list abort) is unchanged - only the namelist placement is corrected.
