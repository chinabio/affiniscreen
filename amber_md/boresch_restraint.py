"""amber_md v2.5.17 — Boresch restraint generation + six-DOF measurement.

PORTED DIRECTLY from the published FEP-SPell-ABFE workflow
(abfe/abfe_alchemy_md.py::BoreschRestraint, abfe/md/cpptraj_tools.py,
abfe/md/analytic.py) so the simulated restraint potential, the measured
reference geometry, and the analytic standard-state correction are GUARANTEED
to use the SAME six internal coordinates.

This replaces amber_md's FEPSetup._write_boresch_RST, whose angle/dihedral atom
orderings did NOT match the coordinates it measured/corrected (see audit).

The six Boresch DOF, with L1/L2/L3 = three ligand atoms and P1/P2/P3 = three
receptor atoms (Boresch 2003; Deng & Roux JCTC 2006):

    r      : distance  L1-P1
    alpha  : angle     P1-L1-L2
    theta  : angle     P2-P1-L1
    gamma  : dihedral  P1-L1-L2-L3
    beta   : dihedral  P2-P1-L1-L2
    phi    : dihedral  P3-P2-P1-L1
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
# Portions ported/adapted from FEP-SPell-ABFE (freeenergylab, MIT License).
from __future__ import annotations
import os
import subprocess
from pathlib import Path

# --- cpptraj template (verbatim from FEP-SPell cpptraj_tools.SIX_DOFS_INPUT) ---
SIX_DOFS_INPUT = """parm PRMTOP
trajin MDCRD

distance bnd_r L1 P1 out bnd_r.dat
angle bnd_alpha P1 L1 L2 out bnd_alpha.dat
angle bnd_theta P2 P1 L1 out bnd_theta.dat
dihedral bnd_gamma P1 L1 L2 L3 out bnd_gamma.dat
dihedral bnd_beta P2 P1 L1 L2 out bnd_beta.dat
dihedral bnd_phi P3 P2 P1 L1 out bnd_phi.dat
"""


# --------------------------------------------------------------------------
# 1. MEASURE the six DOF from an equilibrated trajectory (FEP-SPell logic)
# --------------------------------------------------------------------------
def _mean_std(datfile, is_periodic):
    """Mean/std of a cpptraj .dat column. Circular mean for periodic dihedrals
    (matches cpptraj_tools.mean_std)."""
    import math
    vals = []
    for line in Path(datfile).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        vals.append(float(line.split()[1]))
    n = len(vals)
    if n == 0:
        raise ValueError(f"no data in {datfile}")
    if is_periodic:
        s = sum(math.sin(math.radians(v)) for v in vals) / n
        c = sum(math.cos(math.radians(v)) for v in vals) / n
        mean = math.degrees(math.atan2(s, c))
        R = math.hypot(s, c)
        std = math.degrees(math.sqrt(max(0.0, -2.0 * math.log(R)))) if R > 0 else 0.0
    else:
        mean = sum(vals) / n
        std = (sum((v - mean) ** 2 for v in vals) / n) ** 0.5
    return mean, std


def measure_six_dofs(prmtop, traj, lig_atoms, rec_atoms, workdir,
                     amberbin="", cpptraj=None):
    """Run cpptraj to measure r, alpha, theta, gamma, beta, phi from `traj`.

    lig_atoms = [L1, L2, L3], rec_atoms = [P1, P2, P3] as Amber masks
    (e.g. ":1@C4"). Returns dict of means + the per-DOF std for QC.
    """
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    L1, L2, L3 = lig_atoms
    P1, P2, P3 = rec_atoms
    traj_in = (SIX_DOFS_INPUT
               .replace("PRMTOP", str(prmtop)).replace("MDCRD", str(traj))
               .replace("L1", L1).replace("L2", L2).replace("L3", L3)
               .replace("P1", P1).replace("P2", P2).replace("P3", P3))
    (workdir / "traj.in").write_text(traj_in)
    exe = cpptraj or os.path.join(amberbin, "cpptraj") if amberbin else (cpptraj or "cpptraj")
    subprocess.run([exe, "-i", "traj.in"], cwd=str(workdir), check=True)
    r0, r0s        = _mean_std(workdir / "bnd_r.dat",     False)
    alpha0, a_s    = _mean_std(workdir / "bnd_alpha.dat", False)
    theta0, t_s    = _mean_std(workdir / "bnd_theta.dat", False)
    gamma0, g_s    = _mean_std(workdir / "bnd_gamma.dat", True)
    beta0, b_s     = _mean_std(workdir / "bnd_beta.dat",  True)
    phi0, p_s      = _mean_std(workdir / "bnd_phi.dat",   True)
    refs = dict(r0=r0, alpha0=alpha0, theta0=theta0,
                gamma0=gamma0, beta0=beta0, phi0=phi0)
    stds = dict(r0=r0s, alpha0=a_s, theta0=t_s, gamma0=g_s, beta0=b_s, phi0=p_s)
    (workdir / "boresch_six_dofs_stats_info.dat").write_text(
        "      {:>9}{:>9}{:>9}{:>9}{:>9}{:>9}\n".format(
            "r", "alpha", "theta", "gamma", "beta", "phi") +
        "mean: {r0:9.3f}{alpha0:9.3f}{theta0:9.3f}{gamma0:9.3f}{beta0:9.3f}{phi0:9.3f}\n".format(**refs) +
        " std: {:9.3f}{:9.3f}{:9.3f}{:9.3f}{:9.3f}{:9.3f}\n".format(
            r0s, a_s, t_s, g_s, b_s, p_s))
    return refs, stds


# --------------------------------------------------------------------------
# 2. RESOLVE masks -> 1-based prmtop atom indices (parmed; FEP-SPell logic)
# --------------------------------------------------------------------------
def resolve_atom_indices(prmtop, lig_masks, rec_masks, lig_idx=0):
    """[L1,L2,L3],[P1,P2,P3] masks -> ([l1,l2,l3],[p1,p2,p3]) 1-based indices."""
    import parmed as pmd
    mol = pmd.load_file(str(prmtop))
    lig = []
    for m in lig_masks:
        sel = [i + 1 for i in pmd.amber.AmberMask(mol, m).Selected()]
        lig.append(sel[lig_idx])
    rec = []
    for m in rec_masks:
        sel = [i + 1 for i in pmd.amber.AmberMask(mol, m).Selected()]
        rec.append(sel[0])
    return lig, rec


# --------------------------------------------------------------------------
# 3. WRITE restraints.inp / boresch.RST  (verbatim FEP-SPell atom ordering)
# --------------------------------------------------------------------------
def gen_restraint_str(lig_idx_atoms, rec_idx_atoms, refs, k_bond, k_angle, k_dih):
    """Build the six &rst records EXACTLY as FEP-SPell BoreschRestraint._gen_rst_str.

    lig_idx_atoms = [L1,L2,L3], rec_idx_atoms = [P1,P2,P3]  (1-based ints)
    refs = dict(r0, alpha0, theta0, gamma0, beta0, phi0)  (degrees / angstrom)
    """
    dis_fmt = "&rst iat=%d,%d r1=0.0, r2=%.2f, r3=%.2f, r4=99.0, rk2=%.2f, rk3=%.2f/"
    ang_fmt = "&rst iat=%d,%d,%d r1=0.0, r2=%.2f, r3=%.2f, r4=180.0, rk2=%.2f, rk3=%.2f/"
    dih_fmt = "&rst iat=%d,%d,%d,%d r1=%.2f, r2=%.2f, r3=%.2f, r4=%.2f, rk2=%.2f, rk3=%.2f/"
    L1, L2, L3 = lig_idx_atoms
    P1, P2, P3 = rec_idx_atoms
    r0 = refs["r0"]; alpha0 = refs["alpha0"]; theta0 = refs["theta0"]
    gamma0 = refs["gamma0"]; beta0 = refs["beta0"]; phi0 = refs["phi0"]
    KR = k_bond; KALPHA = KTHETA = k_angle; KGAMMA = KBETA = KPHI = k_dih
    dih11, dih14 = gamma0 - 180.0, gamma0 + 180.0
    dih21, dih24 = beta0 - 180.0, beta0 + 180.0
    dih31, dih34 = phi0 - 180.0, phi0 + 180.0
    rst = [
        dis_fmt % (L1, P1, r0, r0, KR, KR),
        ang_fmt % (P1, L1, L2, alpha0, alpha0, KALPHA, KALPHA),
        ang_fmt % (P2, P1, L1, theta0, theta0, KTHETA, KTHETA),
        dih_fmt % (P1, L1, L2, L3, dih11, gamma0, gamma0, dih14, KGAMMA, KGAMMA),
        dih_fmt % (P2, P1, L1, L2, dih21, beta0, beta0, dih24, KBETA, KBETA),
        dih_fmt % (P3, P2, P1, L1, dih31, phi0, phi0, dih34, KPHI, KPHI),
    ]
    return "\n".join(rst)


def write_restraints_inp(path, lig_idx_atoms, rec_idx_atoms, refs,
                         k_bond=10.0, k_angle=100.0, k_dih=100.0):
    txt = gen_restraint_str(lig_idx_atoms, rec_idx_atoms, refs,
                            k_bond, k_angle, k_dih)
    Path(path).write_text(txt + "\n")
    return txt


# --------------------------------------------------------------------------
# 4. ANALYTIC correction (Deng & Roux Eq.38/40; verbatim FEP-SPell analytic.py)
# --------------------------------------------------------------------------
def restraint_correction(T, r0, theta0, kr, ktheta, kalpha, kphi, kbeta, kgamma):
    from math import pi, log, sin, sqrt
    T = float(T); r0 = float(r0); theta0 = float(theta0)
    kr *= 2.0; ktheta *= 2.0; kalpha *= 2.0; kphi *= 2.0; kbeta *= 2.0; kgamma *= 2.0
    V = 1660.0
    kB = 0.0019872041
    FtC0 = -(kB * T) * log(r0 ** 2 * sin(theta0 / 180.0 * pi)
                           * sqrt((2 * pi * kB * T) ** 3 / (kr * ktheta * kphi)) / V)
    Fr = -(kB * T) * log(1.0 / (8 * pi ** 2)
                         * sqrt((2 * pi * kB * T) ** 3 / (kalpha * kbeta * kgamma)))
    return FtC0 + Fr


def correction_from_refs(refs, T, k_bond=10.0, k_angle=100.0, k_dih=100.0):
    """Convenience: analytic dG_release from measured refs + the SAME force
    constants used in the restraint file."""
    return restraint_correction(
        T=T, r0=refs["r0"], theta0=refs["theta0"],
        kr=k_bond, ktheta=k_angle, kalpha=k_angle,
        kphi=k_dih, kbeta=k_dih, kgamma=k_dih)


if __name__ == "__main__":
    # Self-test 1: analytic correction reproduces the published 11.62 kcal/mol.
    val = restraint_correction(T=298.15, r0=4.470, theta0=101.230,
                               kr=10, ktheta=100, kalpha=100,
                               kphi=100, kbeta=100, kgamma=100)
    print(f"[1] analytic dG_release = {val:.3f} kcal/mol (ref 11.62)")
    assert abs(val - 11.62) < 0.05, "analytic correction drift!"

    # Self-test 2: restraint file uses the FEP-SPell atom ordering exactly.
    refs = dict(r0=4.47, alpha0=72.14, theta0=101.23,
                gamma0=77.887, beta0=-137.094, phi0=9.048)
    txt = gen_restraint_str([11, 12, 13], [101, 102, 103], refs, 10, 100, 100)
    print("[2] restraints.inp:")
    print(txt)
    # the angle/dihedral records must start with P (receptor) atoms, per FEP-SPell
    lines = txt.splitlines()
    assert lines[1].startswith("&rst iat=101,11,12"), "alpha ordering wrong"
    assert lines[2].startswith("&rst iat=102,101,11"), "theta ordering wrong"
    assert lines[3].startswith("&rst iat=101,11,12,13"), "gamma ordering wrong"
    assert lines[4].startswith("&rst iat=102,101,11,12"), "beta ordering wrong"
    assert lines[5].startswith("&rst iat=103,102,101,11"), "phi ordering wrong"
    print("\nALL SELF-TESTS PASSED")
