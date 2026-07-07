# v2.5.0 build final44: sane Amber MM-GBSA default + GUI-driven ABFE per-window time

## 1. Amber MM-GBSA production default 50 ns -> 10 ns
50 ns is excessive for MM-GBSA (an endpoint method; GB energy averages converge
quickly). Lowered to a sane 10 ns in BOTH places so they agree:
  * run_amber.py:  --prod-ns default 50.0 -> 10.0
  * GUI MM-GBSA "Production MD (ns)" Amber default 50 -> 10
OpenMM MM-GBSA stays 5 ns. Equilibration stays 1 ns.
(CLI users who explicitly pass --prod-ns are unaffected.)

## 2. ABFE/RBFE per-window MD time now driven by the GUI
Previously the GUI collected complex_ns/solvent_ns but the Amber ABFE branch
NEVER passed them -> per-window sampling silently used the driver default
(nstlim_prod=1,000,000 = 2 ns/window). Now the Amber ABFE command converts the
field to driver flags:
    nstlim = ns * 1e6 / 2            (dt = 0.002 ps)
    complex_ns      -> --nstlim-prod     (per-window production; both legs share
                                          the driver's single per-window value)
    abfe_equil_ns   -> --nstlim-eq       (default 0.5 ns if unset)
Falls back to the driver default when complex_ns is absent, so behaviour is
unchanged unless the user sets the field.

Note: the Amber fep_driver exposes ONE per-window --nstlim-prod, so complex and
solvent legs use the same per-window length; solvent_ns remains informational
for the compute estimate. (OpenMM/OpenFE RBFE/ABFE already consumed
complex_ns via sim_time_ns -- unchanged.)

## Conversion sanity
   2 ns -> 1,000,000 steps   (== old config default, confirms dt=2 fs)
   5 ns -> 2,500,000
  10 ns -> 5,000,000
  0.5ns ->   250,000

## Scope / safety
GUI + run_amber.py default only. No change to the running ABFE production job
(launched via CLI with explicit --nstlim-prod 2,500,000) or to the analysis
fixes (final38-43). Cumulative.

## Verified
* run_amber.py + GUI compile.
* ns->nstlim conversion validated.
