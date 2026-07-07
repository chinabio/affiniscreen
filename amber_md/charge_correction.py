"""Poisson-Boltzmann finite-size charge correction for ABFE of charged ligands.

Adapted for the Amber MD / OpenFE workflow from the FEP-SPell-ABFE reference
implementation (freeenergylab, MIT License) and the underlying literature:

  (1) Rocklin, Mobley, Dill, Hunenberger, "Calculating the Binding Free
      Energies of Charged Species ... An Accurate Correction Scheme for
      Electrostatic Finite-Size Effects", J. Chem. Phys. 139, 184103 (2013).
  (2) Chen, Deng, Russell, Wu, Abel, Wang, "Accurate Calculation of Relative
      Binding Free Energies between Ligands with Different Net Charges",
      J. Chem. Theory Comput. 14, 6346 (2018).

WHY THIS MATTERS FOR ABFE
-------------------------
When a ligand carries a net charge, the alchemical decoupling/decharging leg
under periodic boundary conditions with lattice-sum (PME) electrostatics
suffers a finite-size artifact that does NOT cancel between the complex and
solvent legs (because the box, the surroundings and the dielectric response
differ). Without this correction an absolute binding free energy for a
charged ligand can be off by several kcal/mol. The v2.5.15 amber_md workflow
had NO such correction; this module supplies it.

The correction has three analytic components (Rocklin 2013):
  NET_USV : net-charge undersolvation in a finite box       (Eq. 15 + 16)
  RIP     : residual integrated potential                   (Eq. 17, via APBS)
  DSC     : discrete-solvent / quadrupole-trace term        (Eq. 35)
The empirical (EMP) term is negligible (Table IV of Ref. 1) and omitted.

REQUIREMENTS
------------
  * AmberTools (cpptraj) on PATH or via $AMBERHOME/bin
  * APBS on PATH                (https://apbs.readthedocs.io)
  * parmed, numpy, pandas, gridDataFormats (`pip install gridDataFormats`)

This module is import-safe: heavy/optional deps are imported lazily inside the
functions so the rest of the workflow can import it without APBS installed.
"""

from __future__ import annotations
import math
import os
import sys
import time

# --- AMBER bin discovery (works with or without $AMBERHOME) ---------------
def _amber_bin(exe):
    home = os.environ.get("AMBERHOME")
    if home:
        cand = os.path.join(home, "bin", exe)
        if os.path.exists(cand):
            return cand
    return exe  # rely on PATH


def align_complex(prmtop, mdcrd, solvent_mask, output_file="align.cpptraj.in"):
    """Autoimage + center + principal-axis align the last frame for APBS."""
    solute_mask = f"!{solvent_mask}"
    basename, ext = os.path.splitext(os.path.basename(mdcrd))
    new_mdcrd = f"{basename}.autoimage{ext}"
    with open(output_file, "w") as fo:
        fo.write("parm %s\n" % prmtop)
        fo.write("trajin %s lastframe\n" % mdcrd)
        fo.write("autoimage anchor %s\n" % solute_mask)
        fo.write("center %s origin\n" % solute_mask)
        fo.write("principal '%s & !@H=' dorotation\n" % solute_mask)
        fo.write("trajout %s\n" % new_mdcrd)
        fo.write("go\nquit\n")
    status = os.system(" ".join([_amber_bin("cpptraj"), "-i", output_file]))
    if status != 0:
        print(time.strftime("%c"))
        print("cpptraj failed to align system during charge correction!")
        sys.exit(1)
    return new_mdcrd


def compute_charge_correction(prmtop, mdcrd, lig_rname, solvent_mask,
                              temperature, epsilon_solv, wat_rname="WAT"):
    """Compute the PB finite-size charge correction (kcal/mol).

    Args:
        prmtop (str): Amber parm7 file (solvated complex OR solvated ligand).
        mdcrd (str): aligned trajectory/restart with box info.
        lig_rname (str): ligand residue name (e.g. "MOL" or "LIG").
        solvent_mask (str): Amber mask of solvent + counter-ions to strip,
            e.g. ":WAT,K+,Na+,Cl-".
        temperature (float): simulation temperature (K).
        epsilon_solv (float): solvent relative permittivity (97.0 for TIP3P).
        wat_rname (str): water residue name.

    Returns:
        pandas.DataFrame with NET_USV, RIP, DSC and Total (kcal/mol).
    """
    import pandas as pd

    mol = PBCAmberParm(prmtop, xyz=mdcrd)
    # delete extra copies of the ligand if present
    lig_indices = [r.number + 1 for r in mol.residues if r.name == lig_rname]
    if not lig_indices:
        raise ValueError(f"Ligand resname '{lig_rname}' not found in {prmtop}.")
    if len(lig_indices) > 1:
        mol.strip(":" + ",".join(str(i) for i in lig_indices[1:]))

    nwat = sum(1 for r in mol.residues if r.name == wat_rname)
    gamma = mol.compute_quadrupole_trace(wat_rname)
    tol = 1e-3
    lig_atoms = mol.residues[(lig_indices[0] - 1)].atoms
    lig_charge = round(sum(a.charge for a in lig_atoms) / tol) * tol

    mol.strip(solvent_mask)
    is_complex = (len(mol.residues) > 1)
    if not is_complex:
        rip_P = 0.0
        rip_I = compute_apbs_rip(mol, "!:%s" % lig_rname, lig_charge,
                                 temperature, epsilon_solv)
    else:
        rec_charge = round(sum(a.charge for a in mol.atoms) / tol) * tol - lig_charge
        rip_P = compute_apbs_rip(mol, ":%s" % lig_rname, rec_charge,
                                 temperature, epsilon_solv)
        rip_I = compute_apbs_rip(mol, "!:%s" % lig_rname, lig_charge,
                                 temperature, epsilon_solv)

    volume = mol.volume
    COULOMB = APBS.COULOMB
    xi_LS = APBS.xi_LS
    NET_USV = -(xi_LS * COULOMB * (lig_charge ** 2)
                / (2 * epsilon_solv * volume ** (1 / 3)))   # Eq.15+16
    RIP = ((rip_P + rip_I) * lig_charge) / volume            # Eq.17
    DSC = -2 * math.pi * COULOMB * gamma * nwat * lig_charge / (3. * volume)  # Eq.35
    Total = NET_USV + RIP + DSC

    print("NET_USV: % 8.3f RIP: % 8.3f DSC: % 8.3f ChargeCorrection: % 8.2f"
          % (NET_USV, RIP, DSC, Total))

    return pd.DataFrame({
        "lig_net_charge": [lig_charge],
        "NET_USV": [NET_USV],
        "RIP": [RIP],
        "DSC": [DSC],
        "Total (kcal/mol)": [Total],
    })


class PBCAmberParm:
    """Thin wrapper giving an AmberParm a `.volume` and helper methods.

    Subclassing parmed.amber.AmberParm is done lazily so importing this module
    does not require parmed at import time.
    """
    def __new__(cls, *args, **kwargs):
        import parmed as pmd

        class _Impl(pmd.amber.AmberParm):
            @property
            def volume(self):
                volume = self.box[:3].prod()
                if self.pointers["IFBOX"] == 1:
                    return volume
                elif self.pointers["IFBOX"] == 2:
                    return volume * 0.7698004  # truncated octahedron
                return volume

            def compute_orthorhombic_box(self, iso=False):
                import numpy as np
                if self.pointers["IFBOX"] == 1:
                    return self.box[:3]
                elif self.pointers["IFBOX"] >= 2:
                    if iso:
                        return np.ones(3) * self.volume ** (1 / 3.)
                    dims = self.coordinates.max(axis=0) - self.coordinates.min(axis=0)
                    scale = dims / dims[0]
                    lens = scale * (self.volume / scale.prod()) ** (1 / 3.)
                    while np.any(dims > lens):
                        lens *= 1.05
                    return lens

            def zero_charge(self, mask):
                for i in pmd.amber.AmberMask(self, mask).Selected():
                    self.atoms[i].charge *= 0.0

            def compute_quadrupole_trace(self, resname="WAT"):
                mask = pmd.amber.AmberMask(self, f":{resname}")
                resid = self.atoms[next(mask.Selected())].residue.number + 1
                mask = pmd.amber.AmberMask(self, ":%d" % resid)
                centers = [i for i in mask.Selected() if self.atoms[i].epsilon > 0.0]
                if len(centers) != 1:
                    raise RuntimeError(
                        "Can only compute quadrupole trace for a water model "
                        "with one vdW center.")
                xj = self.coordinates[centers[0]]
                gamma = 0.0
                for i in mask.Selected():
                    xi = self.coordinates[i]
                    gamma += self.atoms[i].charge * ((xi - xj) ** 2).sum()
                return gamma

        return _Impl(*args, **kwargs)


def compute_apbs_rip(mol, chg_mask, net_charge, temperature, epsilon_solv):
    """Compute the RIP energy term using APBS (kcal-A^3/mol-e, AMBER units)."""
    tmp_mol = mol[:]
    tmp_mol.zero_charge(chg_mask)
    tmp_mol.save("apbs.pqr", None, True)
    del tmp_mol
    is_complex = (len(mol.residues) > 1)
    box = mol.compute_orthorhombic_box(iso=(not is_complex))
    apbs = APBS("apbs.pqr", temperature, epsilon_solv, box)
    apbs.run("apbs.in")
    return apbs.compute_rip("apbs.dx", net_charge)


class APBS:
    """Create and launch APBS single-point Poisson-Boltzmann jobs."""
    xi_CB = -2.380077            # cubic box (Ref.1 Sec. II D)
    xi_LS = -2.837297            # cubic box (Ref.1 Sec. II D)
    BOLTZMANN = 0.001987192      # kcal/mol-K
    AMBER_ELECTROSTATIC = 18.2223
    COULOMB = AMBER_ELECTROSTATIC ** 2   # kcal-A/mol

    def __init__(self, pqr, temperature, sdie, box, grid_spacing=0.2, pdie=1.0):
        import numpy as np
        self._pqr = str(pqr)
        self.temperature = max(0.0, float(temperature))
        self.sdie = max(1.0, float(sdie))
        self.box = np.array(box, dtype=np.float64).flatten()
        assert self.box.size == 3
        self.grid_spacing = max(0.01, float(grid_spacing))
        self.pdie = max(1.0, float(pdie))

    @property
    def dime(self):
        import numpy as np
        ngrid = lambda k: 2 ** k + 1
        k = 6
        dime = ngrid(k)
        while not np.all((self.box / dime) <= self.grid_spacing) and k < 8:
            k += 1
            dime = ngrid(k)
        return [dime] * 3

    def run(self, apbs_in, do_force=False):
        lines = [
            "read", "    mol pqr %s" % self._pqr, "end", "",
            "elec name potential", "    mg-manual",
            "    dime %s" % (" ".join("%d" % n for n in self.dime)),
            "    glen %s" % (" ".join("%f" % b for b in self.box)),
            "    gcent 0.0 0.0 0.0", "    mol 1", "    lpbe", "    bcfl mdh",
            "    pdie %f" % self.pdie, "    sdie %f" % self.sdie,
            "    chgm spl4", "    srfm smol", "    srad 1.4", "    swin 0.3",
            "    sdens 40.0", "    temp %f" % self.temperature,
            "    calcenergy total",
            "    calcforce %s" % ("total" if do_force else "no"),
            "    write pot dx apbs", "end",
        ]
        with open(apbs_in, "w") as f:
            f.write("\n".join(lines))
        os.system(" ".join(["apbs", apbs_in, ">", "apbs.out"]))

    def compute_rip(self, dx_file, net_charge):
        import gridData
        dxfile = gridData.Grid(dx_file)
        dV = dxfile.delta.prod()
        V = dxfile.grid.size * dV
        B_X = dxfile.grid.sum() * dV * self.BOLTZMANN * self.temperature   # Eq.19
        B_QX = -(self.xi_CB * self.COULOMB / self.sdie
                 * net_charge * V ** (2. / 3.))                            # Eq.21
        return B_X - B_QX                                                  # Eq.18


# Convenience: relative permittivity of common water models (for sdie)
EPSILON_SOLV = {"tip3p": 97.0, "tip4pew": 63.0, "opc": 78.4, "spce": 68.0}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="PB finite-size charge correction for a charged ligand.")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--mdcrd", required=True)
    parser.add_argument("--lig-rname", default="MOL")
    parser.add_argument("--solvent-mask", default=":WAT,K+,Na+,Cl-")
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--water", default="tip3p",
                        choices=sorted(EPSILON_SOLV))
    parser.add_argument("--out", default="charge_correction.csv")
    args = parser.parse_args()

    work_dir = os.path.dirname(os.path.abspath(args.mdcrd)) or "."
    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        new_mdcrd = align_complex(args.prmtop, args.mdcrd, args.solvent_mask)
        df = compute_charge_correction(
            prmtop=args.prmtop, mdcrd=new_mdcrd, lig_rname=args.lig_rname,
            solvent_mask=args.solvent_mask, temperature=args.temperature,
            epsilon_solv=EPSILON_SOLV[args.water], wat_rname="WAT")
        df.to_csv(args.out, index=False, float_format="%.4f")
        print(f"Wrote {os.path.join(work_dir, args.out)}")
    finally:
        os.chdir(cwd)
