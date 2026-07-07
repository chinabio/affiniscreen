
# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .utils import run, which_or_die, ensure_dir
from .logger import get_logger
from .config import SystemConfig
log = get_logger()
_WATER_MAP = {"tip3p":("leaprc.water.tip3p","TIP3PBOX"),
              "spce":("leaprc.water.spce","SPCBOX"),
              "tip4pew":("leaprc.water.tip4pew","TIP4PEWBOX"),
              "opc":("leaprc.water.opc","OPCBOX")}

class SystemBuilder:
    def __init__(self, work_dir, cfg):
        self.work_dir = ensure_dir(work_dir/"build")
        self.cfg = cfg

    def _solvate_cmd(self, buf):
        _, box = _WATER_MAP[self.cfg.water_model.lower()]
        kw = "solvateOct" if self.cfg.box_shape == "octahedral" else "solvateBox"
        return f"{kw} system {box} {buf}"

    def _ion_cmds(self):
        cmds = []
        method = (self.cfg.ion_method or "rand").lower()
        n_salt = int(round(self.cfg.salt_conc_M * 100)) if self.cfg.salt_conc_M > 0 else 0

        if method == "rand":
            if self.cfg.neutralize and n_salt > 0:
                cmds.append(f"addIonsRand system Na+ {n_salt} Cl- {n_salt}")
            elif self.cfg.neutralize:
                cmds.append("addIonsRand system Na+ 0 Cl- 0")
            elif n_salt > 0:
                cmds.append(f"addIonsRand system Na+ {n_salt} Cl- {n_salt}")
        elif method == "grid":
            if self.cfg.neutralize:
                cmds.append("addIons2 system Na+ 0")
                cmds.append("addIons2 system Cl- 0")
            if n_salt > 0:
                cmds.append(f"addIonsRand system Na+ {n_salt} Cl- {n_salt}")
        else:
            raise ValueError(
                f"Unknown ion_method {method!r}. Use 'rand' (fast) or 'grid' (slow).")
        return cmds

    def build(self, protein_pdb, ligand_mol2, ligand_frcmod, resname="LIG"):
        which_or_die("tleap")
        water_leaprc, _ = _WATER_MAP[self.cfg.water_model.lower()]
        prmtop = self.work_dir/"complex.prmtop"
        inpcrd = self.work_dir/"complex.inpcrd"
        outpdb = self.work_dir/"complex_solv.pdb"
        leap_in = self.work_dir/"tleap.in"

        ion_block = chr(10).join(self._ion_cmds())
        log.info("tleap ion method: %s  =>  %s",
                 self.cfg.ion_method, ion_block.replace(chr(10), " ; "))

        leap_in.write_text(f"""source leaprc.protein.{self.cfg.protein_ff}
source leaprc.{self.cfg.ligand_ff}
source {water_leaprc}
loadAmberParams {ligand_frcmod.resolve()}
{resname} = loadMol2 {ligand_mol2.resolve()}
protein = loadPDB {protein_pdb.resolve()}
system = combine {{ protein {resname} }}
{self._solvate_cmd(self.cfg.box_buffer_A)}
{ion_block}
saveAmberParm system {prmtop} {inpcrd}
savePDB system {outpdb}
quit
""")
        run(["tleap", "-f", str(leap_in)], cwd=self.work_dir)
        if not prmtop.exists():
            raise RuntimeError(
                "tleap failed; see leap.log. Common causes: missing protein "
                "atoms, unrecognised residue names, or charge mismatch.")
        return prmtop, inpcrd, outpdb
