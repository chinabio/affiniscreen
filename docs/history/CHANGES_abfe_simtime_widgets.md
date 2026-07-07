# v2.5.0 build final45: finish ABFE per-window MD-time exposure

## Decisions (made on request)
1. Independent complex/solvent per-window nstlim: DECLINED. The Amber fep_driver
   exposes a single per-window --nstlim-prod; threading two values through the
   leg builder + LSF templates is a deep, risky change for marginal benefit
   while a production run is live. Complex/solvent share one per-window length.
2. Config-level nstlim_prod default: LEFT at 2 ns/window for CLI back-compat
   (no safe universal change). Instead the GUI default is made sane.

## Changes (gui/pages/0_Setup_and_Launch.py)
* final44 wired complex_ns -> --nstlim-prod and read abfe_equil_ns, but no GUI
  widget existed for equilibration. Added it:
    - "Production / window (ns)"     -> complex_ns      (default 5 ns)
    - "Equilibration / window (ns)"  -> abfe_equil_ns   (default 0.5 ns)
* Relabelled to PER-WINDOW (the value now drives --nstlim-prod, not a total),
  with help text spelling out nstlim = ns x 1e6 / 2 and that legs share it.
* complex_ns default 10 -> 5 ns/window (sane production value; with ~38 dense
  windows 10 ns/window is a large bill).
* solvent_ns kept for the compute estimate but mirrors complex_ns.

## Scope / safety
GUI only. No engine/CLI default change; running ABFE job unaffected. Cumulative
over final38-44. Compiles.
