"""Boresch-style restraint selection + analytical RRHO correction.

Reference: Boresch et al., J. Phys. Chem. B 107, 9535 (2003), Eq. 32.
Cross-checked against Deng & Roux, J. Chem. Theory Comput. 2, 1255 (2006)
and the FEP-SPell-ABFE reference implementation (validated value
11.62 kcal/mol for the T4-lysozyme test restraint).

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.
#
# v2.5.16 FIX (correctness): the analytical Boresch correction is now
# physically consistent with the restraint potential actually written to
# the Amber &rst file. Amber NMR restraints use U = rk*(x - x0)^2, so the
# harmonic spring constant is k = 2*rk. The previous implementation plugged
# the raw rk values into the closed form, producing an error of
# sqrt(2^6) = 8 inside the logarithm, i.e. a systematic ~+1.2 kcal/mol bias
# on EVERY absolute binding free energy. The fix doubles all six force
# constants and is validated against the published reference value.
"""

from __future__ import annotations
import math
from pathlib import Path

KB_KCAL = 0.001987204            # kcal/(mol*K)
V0_A3   = 1660.5392              # A^3 per standard state of 1 M


def select_boresch_atoms(receptor_pdb, ligand_resname,
                         cutoff_A=10.0, target_distance_A=5.0):
    """Pick 3 receptor + 3 ligand atoms suitable for Boresch restraints."""
    try:
        import MDAnalysis as mda
    except ImportError as e:
        raise ImportError(
            "select_boresch_atoms() needs MDAnalysis. "
            "Install with `conda install -c conda-forge MDAnalysis`."
        ) from e
    import numpy as np

    u = mda.Universe(str(receptor_pdb))
    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if len(lig) < 3:
        raise ValueError(f"Ligand {ligand_resname} has <3 heavy atoms.")
    rec = u.select_atoms("protein and name CA")
    if len(rec) < 3:
        raise ValueError("Receptor has <3 Calpha atoms.")

    com = lig.center_of_mass()
    dists = ((lig.positions - com) ** 2).sum(axis=1)
    A = lig[int(dists.argmin())]
    dA = ((lig.positions - A.position) ** 2).sum(axis=1)
    B = lig[int(dA.argmax())]
    AB = B.position - A.position
    cross_norms = np.linalg.norm(np.cross(lig.positions - A.position, AB), axis=1)
    A_local = int(np.where(lig.indices == A.index)[0][0])
    B_local = int(np.where(lig.indices == B.index)[0][0])
    cross_norms[A_local] = -1
    cross_norms[B_local] = -1
    C = lig[int(cross_norms.argmax())]

    dRA = ((rec.positions - A.position) ** 2).sum(axis=1) ** 0.5
    within = np.where(dRA < cutoff_A)[0]
    if len(within) < 3:
        within = np.argsort(dRA)[:3]
    pref = within[np.argsort(np.abs(dRA[within] - target_distance_A))]
    a = rec[int(pref[0])]
    a_idx = int(np.where(rec.indices == a.index)[0][0])
    b_idx = min(len(rec) - 1, a_idx + 3)
    if b_idx == a_idx:
        b_idx = min(len(rec) - 1, a_idx + 1)
    c_idx = min(len(rec) - 1, a_idx + 6) if b_idx == a_idx + 3 else max(0, a_idx - 3)
    if c_idx == a_idx or c_idx == b_idx:
        c_idx = max(0, a_idx - 3)
    b = rec[b_idx]
    c = rec[c_idx]

    r0   = float(np.linalg.norm(A.position - a.position))
    thA0 = _angle_deg(b.position, a.position, A.position)
    thB0 = _angle_deg(a.position, A.position, B.position)
    phA0 = _dihedral_deg(c.position, b.position, a.position, A.position)
    phB0 = _dihedral_deg(b.position, a.position, A.position, B.position)
    phC0 = _dihedral_deg(a.position, A.position, B.position, C.position)

    out = {
        "aA": int(a.index) + 1, "bA": int(b.index) + 1, "cA": int(c.index) + 1,
        "A":  int(A.index) + 1, "B":  int(B.index) + 1, "C":  int(C.index) + 1,
        "r0": r0, "thA0": thA0, "thB0": thB0,
        "phA0": phA0, "phB0": phB0, "phC0": phC0,
        # rk values written to the Amber &rst file (U = rk*(x-x0)^2):
        "kr":  10.0,
        "kth": 100.0,
        "kph": 100.0,
    }
    # v2.5.18: canonical FEP-SPell six-DOF labels so the simulated restraint,
    # the geometry, and the analytic correction all use the SAME coordinates.
    out.update(canonical_dofs_from_legacy(out))
    return out


def canonical_dofs_from_legacy(b):
    """Map the legacy select_boresch_atoms dict to canonical FEP-SPell DOF.

    FEP-SPell coordinates (Boresch 2003; Deng & Roux JCTC 2006):
        r=L1-P1, alpha=P1-L1-L2, theta=P2-P1-L1,
        gamma=P1-L1-L2-L3, beta=P2-P1-L1-L2, phi=P3-P2-P1-L1
    with L1=A,L2=B,L3=C (ligand) and P1=aA,P2=bA,P3=cA (receptor).
    """
    return {
        "L1": int(b["A"]),  "L2": int(b["B"]),  "L3": int(b["C"]),
        "P1": int(b["aA"]), "P2": int(b["bA"]), "P3": int(b["cA"]),
        "alpha0": float(b["thB0"]),   # P1-L1-L2
        "theta0": float(b["thA0"]),   # P2-P1-L1
        "gamma0": float(b["phC0"]),   # P1-L1-L2-L3
        "beta0":  float(b["phB0"]),   # P2-P1-L1-L2
        "phi0":   float(b["phA0"]),   # P3-P2-P1-L1
    }


def _angle_deg(p1, p2, p3):
    import numpy as np
    v1 = p1 - p2
    v2 = p3 - p2
    c = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _dihedral_deg(p1, p2, p3, p4):
    import numpy as np
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = float(np.dot(n1, n2))
    y = float(np.dot(m1, n2))
    return math.degrees(math.atan2(y, x))


def boresch_dG_release(b, T=298.0):
    """Standard-state free energy to RELEASE the 6-DOF Boresch restraint.

    Boresch 2003, Eq. 32:

        dG_release = + kT * ln[  8 pi^2 V0 * sqrt(Kr Kthа Kthb Kpha Kphb Kphc)
                                 -------------------------------------------- ]
                               (r0^2 sin(thA) sin(thB) (2 pi kT)^3)

    where Kx are the HARMONIC spring constants (U = K (x-x0)^2). The values
    stored in `b` ("kr","kth","kph") are the Amber &rst rk constants and the
    Amber restraint energy is U = rk (x-x0)^2, so K = 2*rk for every DOF.

    Returns a positive number (kcal/mol): the cost of removing the restraint
    on the fully-decoupled ligand in the gas-phase standard-state volume.
    """
    kT = KB_KCAL * float(T)
    r0  = float(b["r0"])
    thA = math.radians(float(b["thA0"]))
    thB = math.radians(float(b["thB0"]))

    # Amber rk -> harmonic spring constant K = 2*rk for ALL six DOF.
    Kr   = 2.0 * float(b["kr"])
    KthA = 2.0 * float(b["kth"])
    KthB = 2.0 * float(b["kth"])
    KphA = 2.0 * float(b["kph"])
    KphB = 2.0 * float(b["kph"])
    KphC = 2.0 * float(b["kph"])

    num = (8.0 * math.pi ** 2 * V0_A3
           * math.sqrt(Kr * KthA * KthB * KphA * KphB * KphC))
    den = (r0 ** 2 * math.sin(thA) * math.sin(thB)
           * (2.0 * math.pi * kT) ** 3)
    return kT * math.log(num / den)


def boresch_correction(b, T=298.0):
    """Analytical Boresch correction dG_rest (kcal/mol), ADD direction.

    This is the term added to the raw alchemical (decharge+vdw) dG of the
    complex leg to account for releasing the restraint to the 1 M standard
    state. It equals -dG_release.

    Backward compatible: same call signature and sign convention as the
    previous release, but now numerically correct (validated to the
    published 11.62 kcal/mol reference for the canonical test restraint).
    """
    return -boresch_dG_release(b, T=T)


def verify_boresch_restraint(boresch, prmtop, coord_file):
    """Read a snapshot and print actual instantaneous values vs. targets."""
    try:
        import MDAnalysis as mda
    except ImportError as e:
        raise ImportError("verify_boresch_restraint() needs MDAnalysis.") from e
    import numpy as np

    u = mda.Universe(str(prmtop), str(coord_file))

    def get(i):
        return u.atoms[i - 1]

    aA = get(boresch["aA"])
    bA = get(boresch["bA"])
    cA = get(boresch["cA"])
    A  = get(boresch["A"])
    B  = get(boresch["B"])
    C  = get(boresch["C"])

    r   = float(np.linalg.norm(aA.position - A.position))
    thA = _angle_deg(bA.position, aA.position, A.position)
    thB = _angle_deg(aA.position, A.position, B.position)
    phA = _dihedral_deg(cA.position, bA.position, aA.position, A.position)
    phB = _dihedral_deg(bA.position, aA.position, A.position, B.position)
    phC = _dihedral_deg(aA.position, A.position, B.position, C.position)

    rows = [
        ("distance r ",  r,   boresch["r0"],   "A"),
        ("angle thA  ",  thA, boresch["thA0"], "deg"),
        ("angle thB  ",  thB, boresch["thB0"], "deg"),
        ("dihed phA  ",  phA, boresch["phA0"], "deg"),
        ("dihed phB  ",  phB, boresch["phB0"], "deg"),
        ("dihed phC  ",  phC, boresch["phC0"], "deg"),
    ]
    print(f"{'restraint':12s}  {'actual':>10s}  {'target':>10s}  "
          f"{'|delta|':>9s}  unit")
    print("-" * 60)
    bad = []
    for name, val, tgt, unit in rows:
        delta = abs(val - tgt)
        if unit == "deg" and "dihed" in name:
            d = (val - tgt + 180) % 360 - 180
            delta = abs(d)
        flag = ""
        if unit == "deg" and "angle" in name and (val < 20 or val > 160):
            flag = "  <-- WARNING: near singularity"
            bad.append(name)
        if unit == "A" and val < 1.5:
            flag = "  <-- WARNING: atoms overlap"
            bad.append(name)
        print(f"{name:12s}  {val:10.3f}  {tgt:10.3f}  {delta:9.3f}  {unit}{flag}")
    if bad:
        print("\n[!] Re-pick atoms before running FEP.")
    else:
        print("\n[OK] Boresch geometry looks sane.")
    return not bool(bad)


def _self_test():
    """Validate the closed form against the published FEP-SPell reference.

    FEP-SPell-ABFE / Deng-Roux T4-lysozyme test restraint:
      r0 = 4.470 A, theta = 101.230 deg, rk_r = 10, rk_ang = rk_dih = 100
      -> dG_release = 11.62 kcal/mol  (paper)  /  11.63 (this code)
    """
    b = {"r0": 4.470, "thA0": 101.230, "thB0": 101.230,
         "kr": 10.0, "kth": 100.0, "kph": 100.0}
    dGrel = boresch_dG_release(b, T=298.15)
    ok = abs(dGrel - 11.62) < 0.05
    print(f"[self-test] dG_release = {dGrel:.3f} kcal/mol "
          f"(reference 11.62) -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import argparse, json, sys
    p = argparse.ArgumentParser(
        description="Auto-pick Boresch atoms and/or verify a geometry.")
    p.add_argument("--pdb", type=Path, help="Receptor PDB for atom picking.")
    p.add_argument("--ligand-resname", default="LIG")
    p.add_argument("--out-json", type=Path,
                   help="Write picked atoms + correction to JSON.")
    p.add_argument("--temperature", type=float, default=298.0)
    p.add_argument("--verify", action="store_true",
                   help="Read --boresch-json + --prmtop + --coord and "
                        "print actual values.")
    p.add_argument("--boresch-json", type=Path)
    p.add_argument("--prmtop", type=Path)
    p.add_argument("--coord",  type=Path)
    p.add_argument("--self-test", action="store_true",
                   help="Validate the closed form against the published "
                        "reference value and exit.")
    a = p.parse_args()

    if a.self_test:
        sys.exit(0 if _self_test() else 1)

    if a.verify:
        if not (a.boresch_json and a.prmtop and a.coord):
            sys.exit("--verify needs --boresch-json, --prmtop, --coord")
        b = json.loads(a.boresch_json.read_text())
        ok = verify_boresch_restraint(b, a.prmtop, a.coord)
        sys.exit(0 if ok else 1)

    if a.pdb:
        b = select_boresch_atoms(a.pdb, a.ligand_resname)
        b["dG_correction_kcal_mol"] = boresch_correction(b, T=a.temperature)
        b["dG_release_kcal_mol"] = boresch_dG_release(b, T=a.temperature)
        print(json.dumps(b, indent=2))
        if a.out_json:
            a.out_json.write_text(json.dumps(b, indent=2))
            print(f"\nWrote {a.out_json}")
        print(f"\nBoresch correction at T={a.temperature} K: "
              f"{b['dG_correction_kcal_mol']:+.3f} kcal/mol "
              f"(dG_release = {b['dG_release_kcal_mol']:+.3f})")
    else:
        p.print_help()


_PRECHECK_R_MIN, _PRECHECK_R_MAX = 4.0, 13.0
_PRECHECK_ANG_MIN, _PRECHECK_ANG_MAX = 30.0, 150.0

def _precheck_boresch_dict(b):
    """(ok, reason) geometric acceptance check for a boresch dict. v2.5.17."""
    if not isinstance(b, dict):
        return True, "non-dict boresch (skipped)"
    bad = []
    r0 = b.get("r0")
    if r0 is not None and not (_PRECHECK_R_MIN <= r0 <= _PRECHECK_R_MAX):
        bad.append("r0=%.2f A outside [%.1f,%.1f]" % (r0, _PRECHECK_R_MIN, _PRECHECK_R_MAX))
    for key in ("thA0", "thB0"):
        a = b.get(key)
        if a is not None and not (_PRECHECK_ANG_MIN <= a <= _PRECHECK_ANG_MAX):
            bad.append("%s=%.1f deg outside [%.0f,%.0f]" % (key, a, _PRECHECK_ANG_MIN, _PRECHECK_ANG_MAX))
    return (len(bad) == 0), ("; ".join(bad) if bad else "ok")


def restraint_correction_dengroux(T, r0, theta0, kr, ktheta, kalpha,
                                  kphi, kbeta, kgamma):
    """Deng & Roux JCTC 2006 Eq.38/40 (FEP-SPell analytic.py). rk -> 2*rk.
    Returns dG_release (kcal/mol, positive)."""
    from math import pi, log, sin, sqrt
    T = float(T); r0 = float(r0); theta0 = float(theta0)
    kr = float(kr)*2.0; ktheta = float(ktheta)*2.0; kalpha = float(kalpha)*2.0
    kphi = float(kphi)*2.0; kbeta = float(kbeta)*2.0; kgamma = float(kgamma)*2.0
    V = 1660.0
    kB = 0.0019872041
    FtC0 = -(kB*T)*log(r0**2*sin(theta0/180.0*pi)
                       *sqrt((2*pi*kB*T)**3/(kr*ktheta*kphi))/V)
    Fr = -(kB*T)*log(1.0/(8*pi**2)
                     *sqrt((2*pi*kB*T)**3/(kalpha*kbeta*kgamma)))
    return FtC0 + Fr


def boresch_correction_dengroux(b, T=298.0):
    """Analytic ADD-direction correction (= -dG_release) from canonical refs,
    using the SAME force constants that go into the &rst file."""
    theta0 = float(b.get("theta0", b.get("thA0")))
    dG_release = restraint_correction_dengroux(
        T=T, r0=float(b["r0"]), theta0=theta0,
        kr=float(b["kr"]), ktheta=float(b["kth"]), kalpha=float(b["kth"]),
        kphi=float(b["kph"]), kbeta=float(b["kph"]), kgamma=float(b["kph"]))
    return -dG_release
