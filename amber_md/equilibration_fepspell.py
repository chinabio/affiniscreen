"""amber_md v2.5.17 — staged equilibration (ported from FEP-SPell-ABFE).

Reproduces abfe/abfe_equilibration.py + abfe/md/amber_mdin.py exactly:

  min-1   restrained minimization   (ntr=1, restraint_wt, restraintmask)
  min-2   free minimization         (ntr=0)
  heat-n  NVT heat ladder, restrained, with TEMP_0 ramp  (one per temp step)
  press-n NPT press, restrained                          (one per temp step)
  relax   NPT, FREE (ntr=0), production-length relaxation

Default heat ladder: 5 -> 100 -> 200 -> 298.15 K  (heat_temps), i.e. THREE
heat/press cycles. The gentle ladder + restrained->free staging is what keeps
the solvated box and the complex pocket stable before alchemy starts; skipping
it is a common cause of eq blow-ups.

This produces the per-step mdin files and an ordered run list. It does NOT
submit anything (no subprocess) so it is safe to generate and inspect.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
# Portions ported/adapted from FEP-SPell-ABFE (freeenergylab, MIT License).
from __future__ import annotations
from collections import OrderedDict
from pathlib import Path

# --- mdin templates (verbatim from FEP-SPell abfe/md/amber_mdin.py) ---------
MIN_1 = """Minimisation Stage 1
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,

  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",
 /
"""

MIN_2 = """Minimisation Stage 2
 &cntrl
  imin = 1,
  ntmin = 1,
  maxcyc = 10000,
  ncyc = 5000,

  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,
 /
"""

HEAT = """Heat MD
 &cntrl
  imin = 0,
  irest = 0,
  ntx = 1,
  nstlim = 10000,
  dt = DT,

  ntb = 1,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 0,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,

  nmropt = 1,
 /

  &wt
  type = 'TEMP0',
  istep1 = 0, istep2 = 8000,
  value1 = TEMPI, value2 = TEMP0,
  /
  &wt
  type = 'END',
  /
"""

PRESS = """Press MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = 10000,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 1,
  restraint_wt = RESTRAINT_WT,
  restraintmask = "RESTRAINTMASK",

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 0,
 /
"""

RELAX = """Relax MD
 &cntrl
  imin = 0,
  irest = 1,
  ntx = 5,
  nstlim = NSTLIM,
  dt = DT,

  ntb = 2,
  cut = CUT,
  ntr = 0,
  restraint_wt = 0.0,

  ntp = 1,
  pres0 = 1.0,
  taup = 2.0,
  barostat = 2,

  ntc = 2,
  ntf = 2,

  ntt = 3,
  gamma_ln = 2.0,
  ig = -1,
  tempi = TEMPI,
  temp0 = TEMP0,

  ioutfm = 1,
  ntpr = 2500,
  ntwx = 2500,
 /
"""

DEFAULT_HEAT_TEMPS = [5.0, 100.0, 200.0, 298.15]


def build_equilibration(temperature=298.15, cutoff=9.0, timestep=0.001,
                        relax_length_ns=5.0,
                        restraint_wt=10.0,
                        restraintmask="!:WAT,Cl-,K+,Na+ & !@H=",
                        heat_temps=None):
    """Return an OrderedDict {step_name: mdin_text} for the full ladder.

    Mirrors abfe_equilibration.py: heat dt is fixed at 0.002 ps during the
    ladder; relax uses the production timestep. The final heat_temps entry is
    forced to `temperature` (as in _parse_args).
    """
    heat_temps = list(heat_temps or DEFAULT_HEAT_TEMPS)
    heat_temps[-1] = temperature
    nstlim_relax = int(relax_length_ns / timestep * 1000)

    steps = OrderedDict()
    steps["min-1"] = (MIN_1
                      .replace("CUT", str(cutoff))
                      .replace("RESTRAINT_WT", str(restraint_wt))
                      .replace("RESTRAINTMASK", restraintmask))
    steps["min-2"] = MIN_2.replace("CUT", str(cutoff))
    for n, (tempi, temp0) in enumerate(zip(heat_temps[:-1], heat_temps[1:]), 1):
        steps[f"heat-{n}"] = (HEAT
                              .replace("DT", str(0.001))
                              .replace("CUT", str(cutoff))
                              .replace("TEMPI", str(tempi))
                              .replace("TEMP0", str(temp0))
                              .replace("RESTRAINT_WT", str(restraint_wt))
                              .replace("RESTRAINTMASK", restraintmask))
        steps[f"press-{n}"] = (PRESS
                               .replace("DT", str(0.001))
                               .replace("CUT", str(cutoff))
                               .replace("TEMP0", str(temp0))
                               .replace("RESTRAINT_WT", str(restraint_wt))
                               .replace("RESTRAINTMASK", restraintmask))
    steps["relax"] = (RELAX
                      .replace("NSTLIM", str(nstlim_relax))
                      .replace("DT", str(timestep))
                      .replace("CUT", str(cutoff))
                      .replace("TEMPI", str(temperature))
                      .replace("TEMP0", str(temperature)))
    return steps


def write_equilibration(workdir, **kw):
    """Write each step's mdin to workdir/<step>.in; return the ordered names."""
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    steps = build_equilibration(**kw)
    for name, txt in steps.items():
        (workdir / f"{name}.in").write_text(txt)
    (workdir / "RUN_ORDER.txt").write_text("\n".join(steps) + "\n")
    return list(steps)


if __name__ == "__main__":
    steps = build_equilibration()
    print("Equilibration step order (default 5->100->200->298.15 K ladder):")
    for i, name in enumerate(steps, 1):
        print(f"  {i:2d}. {name}")
    print(f"\nTotal steps: {len(steps)}")
    # sanity: min-1 restrained, min-2 free, relax free, heat has TEMP0 ramp
    assert "ntr = 1" in steps["min-1"] and "ntr = 0" in steps["min-2"]
    assert "ntr = 0" in steps["relax"]
    assert "TEMP0" in steps["heat-1"]
    assert list(steps).count("relax") == 1
    nheat = sum(1 for k in steps if k.startswith("heat-"))
    npress = sum(1 for k in steps if k.startswith("press-"))
    print(f"heat cycles={nheat}, press cycles={npress} (expect 3 each)")
    assert nheat == 3 and npress == 3
    print("\nself-test OK")
