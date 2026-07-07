# v2.5.38 - wire reliability gate into the cycle-closer

Date: 2026-06-20

_leg() now returns the restraint-leg result; the cycle-closer adds complex_restraint to _subres and computes `unreliable = legs with dG_reliable is False`. all_ok now requires len(unreliable)==0, so a physically broken restraint leg (dV/dl runaway / BAR-TI blowup) forces dG_bind UNTRUSTED + exit 1. ABFE_RESULT.txt shows the restraint leg status + reasons; ABFE_RESULT.json gains unreliable_legs, complex_restraint, restraint_reliable.
