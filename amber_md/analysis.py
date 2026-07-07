"""MM/GBSA + cpptraj RMSD/RMSF + plotting (v2.4.4 - adds prepare_inputs)."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .utils import run, which_or_die, ensure_dir
from .logger import get_logger
from .config import MMGBSAConfig
from .topology import TopologySplitter
log = get_logger()


class TrajectoryAnalyzer:
    def __init__(self, work_dir, prmtop, traj):
        self.dir = ensure_dir(work_dir/"analysis")
        self.prm = Path(prmtop)
        self.trj = Path(traj)

    def rmsd_rmsf(self, ref=None):
        which_or_die("cpptraj")
        rmsd_dat = self.dir/"rmsd.dat"
        rmsf_dat = self.dir/"rmsf.dat"
        ref_line = f"reference {ref}\n" if ref else ""
        script = self.dir/"analysis.cpptraj"
        script.write_text(f"""parm {self.prm}
trajin {self.trj}
{ref_line}autoimage
rms RMSD_BB first @CA,C,N out {rmsd_dat} mass
atomicfluct RMSF out {rmsf_dat} @CA byres
run
quit
""")
        run(["cpptraj", "-i", str(script)], cwd=self.dir)
        return rmsd_dat, rmsf_dat

    def plot(self, rmsd, rmsf):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            return Path(), Path()
        d = np.loadtxt(rmsd, comments=("#", "@"))
        f = np.loadtxt(rmsf, comments=("#", "@"))
        p1 = self.dir/"rmsd.png"
        p2 = self.dir/"rmsf.png"
        plt.figure(); plt.plot(d[:, 0], d[:, 1])
        plt.xlabel("Frame"); plt.ylabel("RMSD (Å)"); plt.tight_layout()
        plt.savefig(p1, dpi=140); plt.close()
        plt.figure(); plt.plot(f[:, 0], f[:, 1])
        plt.xlabel("Residue"); plt.ylabel("RMSF (Å)"); plt.tight_layout()
        plt.savefig(p2, dpi=140); plt.close()
        return p1, p2


class MMGBSAAnalyzer:
    def __init__(self, work_dir, cfg, ligand_resname="LIG"):
        self.dir = ensure_dir(work_dir/"mmgbsa")
        self.work_dir = Path(work_dir)
        self.cfg = cfg
        self.ligand_resname = ligand_resname

    # ---------------- input file ----------------
    def _write_input(self):
        f = self.dir/"mmgbsa.in"
        decomp_block = ""
        if self.cfg.decomposition and self.cfg.decomp_residues:
            decomp_block = f"""
&decomp
  idecomp=2, dec_verbose=1, print_res="{self.cfg.decomp_residues}"
/
"""
        f.write_text(f"""MM/GBSA
&general
  startframe={self.cfg.start_frame}, endframe={self.cfg.end_frame or 9999999},
  interval={self.cfg.stride}, verbose=2, keep_files=0
/
&gb
  igb={self.cfg.igb}, saltcon={self.cfg.salt_conc}
/
{decomp_block}
""")
        return f

    def prepare_inputs(self):
        """v2.4.4: Pre-stage everything the in-LSF MM/GBSA block needs.

        - Creates <work_dir>/mmgbsa/.
        - Writes mmgbsa.in (using the cfg values, same as run_serial()).
        - Returns the path to mmgbsa.in.

        Idempotent: safe to call multiple times. Overwrites mmgbsa.in to
        keep it in sync with the current cfg; if you want a hand-tuned
        file preserved, save it under a different name.
        """
        return self._write_input()

    # ---------------- topologies ----------------
    def ensure_topologies(self, solvated_prm, complex_prm=None,
                          receptor_prm=None, ligand_prm=None):
        # v2.2.9: explicit paths take precedence if all exist
        if all(p and Path(p).exists()
               for p in (complex_prm, receptor_prm, ligand_prm)):
            return {"solvated": solvated_prm, "complex": complex_prm,
                    "receptor": receptor_prm, "ligand": ligand_prm}
        # v2.2.9: auto-discover existing topologies in work_dir/topo/
        topo_dir = Path(self.work_dir)/"topo"
        candidates = {"solvated": topo_dir/"solvated.prmtop",
                      "complex":  topo_dir/"complex.prmtop",
                      "receptor": topo_dir/"receptor.prmtop",
                      "ligand":   topo_dir/"ligand.prmtop"}
        if all(p.exists() and p.stat().st_size > 0 for p in candidates.values()):
            log.info("MMGBSAAnalyzer: using existing topologies in %s", topo_dir)
            return candidates
        # Otherwise, run the splitter (idempotent in v2.2.9)
        splitter = TopologySplitter(self.work_dir, self.ligand_resname)
        topos = splitter.split(solvated_prm)
        splitter.sanity_check(topos)
        return topos

    # ---------------- local serial run (fallback) ----------------
    def run_serial(self, solvated_prm, traj,
                   complex_prm=None, receptor_prm=None, ligand_prm=None):
        which_or_die("MMPBSA.py")
        topos = self.ensure_topologies(
            solvated_prm, complex_prm, receptor_prm, ligand_prm)
        inp = self._write_input()
        out = self.dir/"FINAL_RESULTS_MMPBSA.dat"
        cmd = ["MMPBSA.py", "-O", "-i", str(inp), "-o", str(out),
               "-sp", str(topos["solvated"]), "-cp", str(topos["complex"]),
               "-rp", str(topos["receptor"]), "-lp", str(topos["ligand"]),
               "-y", str(traj)]
        if self.cfg.decomposition:
            cmd += ["-do", str(self.dir/"FINAL_DECOMP_MMPBSA.dat")]
        run(cmd, cwd=self.dir)
        return out

    @staticmethod
    def parse_delta_total(results_dat):
        in_delta = False
        for line in Path(results_dat).read_text().splitlines():
            if "Differences (Complex - Receptor - Ligand)" in line:
                in_delta = True
                continue
            if in_delta and line.strip().startswith("DELTA TOTAL"):
                try:
                    return float(line.split()[2])
                except (IndexError, ValueError):
                    return None
        return None