# v2.5.0 build final43: expose MM-GBSA MD simulation time in the GUI

## Question answered
How much MM-GBSA MD time was configured (and hidden)?
  * Amber engine  (run_amber.py -> pmemd): prod = 50 ns  (--prod-ns default 50.0
                   -> prod_nsteps=25,000,000 @ dt=2 fs), equil = 1 ns.
  * OpenMM/OpenFE  (mmgbsa_openmm.py):      prod = 5 ns  (prod_ns=5.0), equil = 1 ns.
Both were HARD defaults; the GUI never set --prod-ns/--equil-ns for either path,
so users could not change MM-GBSA sampling length from the wizard.

## Change (gui/pages/0_Setup_and_Launch.py)
* MM-GBSA settings expander now has two fields:
    - "Production MD (ns)"     -> params["prod_ns"]
    - "Equilibration MD (ns)"  -> params["equil_ns"]
  The Production default adapts to the engine (Amber 50 ns, OpenMM 5 ns), so
  existing behaviour is unchanged unless the user edits the value.
* _build_commands wiring:
    - OpenMM branch:  --prod-ns  P.get("prod_ns", P.get("complex_ns", 5.0))
    - Amber  branch:  now passes --prod-ns / --equil-ns (previously omitted ->
                      always 50 ns). Reads P.get("prod_ns", 50.0)/("equil_ns",1).
* "Review & launch" compute estimate (total_ns) now counts MM-GBSA prod+equil
  (complex_ns/solvent_ns only exist for ABFE/RBFE).

## Scope / safety
GUI only. No change to CLI defaults, the running ABFE production job, or the
analysis fixes (final38-42). Cumulative over final42 (includes the params
NameError fix).

## Verified
* file compiles; both engine branches pass the flags; defaults preserve prior
  behaviour (Amber 50 ns / OpenMM 5 ns) when untouched.
