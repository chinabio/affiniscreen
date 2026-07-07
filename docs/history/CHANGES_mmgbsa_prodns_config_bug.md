# v2.5.3 build final52: FIX — GUI MM-GBSA prod time ignored on the --config path

## Symptom (user-reported, confirmed from a live prod.out)
GUI showed "Production MD = 10 ns" but the Amber MM-GBSA job ran with
nstlim=25,000,000 (= 50 ns). Effective rate ~13 ns/day -> ~90 h job instead of
~18 h.

## Root cause (a real wiring bug, not user error)
1. amber_md/config.py: MDConfig.prod_nsteps default = 25_000_000 (50 ns) -- the
   last surviving copy of the old 50 ns default (final51 fixed the CLI/GUI/batch
   sites but MISSED the dataclass default).
2. gui/amber_config.py build_config(): set cfg.mmgbsa / cfg.fep / protonation
   but NEVER set cfg.md, so the GUI's prod_ns/equil_ns in P were dropped. The
   serialized wizard_config.json therefore carried the 50 ns default.
3. run_amber.py: `if a.config: cfg = WorkflowConfig.load(config)` builds cfg
   ENTIRELY from the file; the CLI --prod-ns (line 209) is only used in the
   `else` branch. The GUI always passes --config, so --prod-ns 10 was ignored.
   Net: 50 ns default won every time via the GUI.

## Fix (three coordinated changes)
* config.py: MDConfig.prod_nsteps 25_000_000 -> 5_000_000 (10 ns default).
* gui/amber_config.py build_config(): map P["prod_ns"]/P["equil_ns"]
  (fallback complex_ns) into cfg.md.prod_nsteps/equil_nsteps
  (nsteps = ns * 1e6 / 2). Now the GUI field is serialized into the config.
* run_amber.py --config branch: if --prod-ns / --equil-ns are passed with a
  NON-default value, they OVERRIDE the loaded config (belt-and-braces; argparse
  defaults 10/1 do not override).

## Verified END-TO-END (imported the package, not just compiled)
* MDConfig() default -> 5,000,000 steps = 10 ns.
* build_config(prod_ns=10) -> cfg.md.prod_nsteps = 5,000,000 (10 ns).
* save()+load() round-trip (simulates --config) preserves 10 ns.
* build_config(prod_ns=25) -> 12,500,000 (25 ns) — arbitrary values flow.

## Action for the user's RUNNING job
The in-flight job is 50 ns (5.3% done, ~3.5 days left). Recommend bkill +
relaunch via GUI with final52 (will now be a true 10 ns job, ~18 h).

## Carried forward
final46 driver hardening, v2.5.3 submit activate fix, final47 anti-storm guard,
final49 report-kit fixes, final50 unified report button, final51 prod-ns
unification (CLI/GUI/batch).
