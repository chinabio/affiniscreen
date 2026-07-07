
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .config import MDConfig
from .utils import ensure_dir

class MDInputWriter:
    def __init__(self, out_dir, cfg): self.out = ensure_dir(Path(out_dir)/"mdin"); self.c = cfg
    def write_min(self):
        f1 = self.out/"min1.in"; f2 = self.out/"min2.in"
        f1.write_text(f"""Min1
&cntrl
 imin=1,maxcyc={self.c.min_maxcyc},ncyc={self.c.min_ncyc},
 ntb=1,ntr=1,cut={self.c.cutoff_A},restraint_wt={self.c.min_restraint_wt},
 restraintmask='!:WAT,Na+,Cl- & !@H='
/
""")
        f2.write_text(f"""Min2
&cntrl
 imin=1,maxcyc={self.c.min_maxcyc},ncyc={self.c.min_ncyc},
 ntb=1,ntr=0,cut={self.c.cutoff_A}
/
""")
        return f1,f2
    def write_heat(self):
        f = self.out/"heat.in"
        f.write_text(f"""Heat
&cntrl
 imin=0,irest=0,ntx=1,nstlim={self.c.heat_nsteps},dt=0.001,
 ntc={self.c.ntc},ntf={self.c.ntf},ntt=3,gamma_ln=2.0,ig=-1,
 tempi={self.c.heat_T_start},temp0={self.c.heat_T_end},
 ntb=1,ntp=0,cut={self.c.cutoff_A},
 ntpr=1000,ntwx=1000,ntwr=10000,
 ntr=1,restraint_wt=5.0,restraintmask='!:WAT,Na+,Cl- & !@H='
/
""")
        return f
    def write_equil(self):
        f = self.out/"equil.in"
        f.write_text(f"""Equil
&cntrl
 imin=0,irest=1,ntx=5,nstlim={self.c.equil_nsteps},dt=0.001,
 ntc={self.c.ntc},ntf={self.c.ntf},ntt=3,gamma_ln=2.0,ig=-1,
 temp0={self.c.temperature_K},ntp=1,pres0={self.c.pressure_bar},taup=2.0,
 ntb=2,cut={self.c.cutoff_A},
 ntpr=1000,ntwx=5000,ntwr=10000,
 ntr=1,restraint_wt={self.c.equil_restraint_wt},restraintmask='@CA,C,N'
/
""")
        return f
    def write_prod(self):
        f = self.out/"prod.in"
        f.write_text(f"""Prod
&cntrl
 imin=0,irest=1,ntx=5,nstlim={self.c.prod_nsteps},dt={self.c.prod_dt_ps},
 ntc={self.c.ntc},ntf={self.c.ntf},ntt=3,gamma_ln=2.0,ig=-1,
 temp0={self.c.temperature_K},ntp=1,pres0={self.c.pressure_bar},taup=2.0,
 ntb=2,cut={self.c.cutoff_A},
 ntpr={self.c.prod_print_freq},ntwx={self.c.prod_print_freq},
 ntwr={self.c.prod_print_freq*10},iwrap=1
/
""")
        return f
    def write_all(self):
        m1,m2 = self.write_min()
        return {"min1":m1,"min2":m2,"heat":self.write_heat(),
                "equil":self.write_equil(),"prod":self.write_prod()}
