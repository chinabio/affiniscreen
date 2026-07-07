"""amber_md v2.5.17 — automatic Boresch six-atom selection.

Implements the auto-selection that FEP-SPell-ABFE leaves as an unimplemented
TODO (its `else` branch writes a broken 'pass' restraint file). Picks three
ligand atoms (L1,L2,L3) and three receptor Calpha atoms (P1,P2,P3) that define
a STIFF, WELL-CONDITIONED Boresch restraint, following the established
heuristics (BAT.py / Yank / Baumann-Gilson) but RMSF-aware:

  * use the equilibration trajectory to rank atoms by fluctuation (low RMSF =
    rigid = good anchor) -- strictly better than a single-frame pick;
  * ligand: heavy atoms, low RMSF, L1 nearest the ligand centre of mass, then
    L2/L3 chosen for separation + non-collinearity;
  * receptor: Calpha atoms in a 5-11 A shell around L1, low RMSF, mutually
    non-collinear and spread out;
  * reject any choice whose angles fall outside [SANITY_ANG_MIN,
    SANITY_ANG_MAX] (the Jacobian sin(theta) term blows up near 0/180 deg).

Manual masks (if provided by the user) ALWAYS override this. If no valid set is
found the selector RAISES (never writes a 'pass' file).

Selection is transparent: writes boresch_atom_selection.dat with the chosen
atoms, their RMSF, the resulting geometry, and pass/fail of every check.

Runtime deps: parmed, numpy, and cpptraj (for RMSF). cpptraj is invoked by the
caller-provided `rmsf_fn` so this module stays subprocess-free and testable.
"""
from __future__ import annotations
import math
from pathlib import Path

# geometric sanity thresholds (degrees / angstrom)
SANITY_ANG_MIN = 30.0
SANITY_ANG_MAX = 150.0
REC_SHELL_MIN = 5.0
REC_SHELL_MAX = 11.0
LIG_MIN_SEP = 1.5      # min pairwise distance between chosen ligand atoms
REC_MIN_SEP = 4.0      # min pairwise distance between chosen receptor atoms
COLLINEAR_EPS = 1e-3


# --------------------------------------------------------------------------
# vector helpers
# --------------------------------------------------------------------------
def _v(a, b):  # b - a
    return (b[0]-a[0], b[1]-a[1], b[2]-a[2])

def _norm(u):
    return math.sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])

def _dist(a, b):
    return _norm(_v(a, b))

def _angle(a, b, c):
    """angle a-b-c in degrees"""
    u = _v(b, a); w = _v(b, c)
    nu, nw = _norm(u), _norm(w)
    if nu == 0 or nw == 0:
        return 0.0
    cosv = max(-1.0, min(1.0, (u[0]*w[0]+u[1]*w[1]+u[2]*w[2])/(nu*nw)))
    return math.degrees(math.acos(cosv))

def _collinear(a, b, c):
    """True if a,b,c are ~collinear (cross product ~0)."""
    u = _v(a, b); w = _v(a, c)
    cx = (u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0])
    return _norm(cx) < COLLINEAR_EPS * max(1.0, _norm(u)*_norm(w))


# --------------------------------------------------------------------------
# main selector
# --------------------------------------------------------------------------
def autoselect_boresch_atoms(prmtop, inpcrd, ligand_resmask=":LIG",
                             receptor_mask="@CA", rmsf=None,
                             report_path=None):
    """Return ((L1,L2,L3),(P1,P2,P3)) as Amber masks like ':1@C4'.

    Parameters
    ----------
    prmtop, inpcrd : str
        Topology + a representative (equilibrated) coordinate set.
    ligand_resmask : str
        Amber residue mask selecting the ligand (default ':LIG').
    receptor_mask : str
        Atom mask for receptor anchor candidates (default '@CA', i.e. Calpha).
    rmsf : dict[int,float] | None
        Optional {atom_index_0based: rmsf_angstrom} from the equil trajectory.
        When None, selection is geometry-only (single frame). Supplying it
        (recommended) makes the pick fluctuation-aware.
    report_path : str | None
        If set, write a human-readable selection report there.
    """
    import parmed as pmd
    import numpy as np

    mol = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    coords = mol.coordinates  # (natom,3)

    def amber_idx(atom):  # 1-based residue + @name mask
        return f":{atom.residue.idx+1}@{atom.name}"

    def rmsf_of(i):
        if rmsf is None:
            return 0.0
        return rmsf.get(i, 0.0)

    # ---- candidate pools ----
    lig_sel = set(pmd.amber.AmberMask(mol, ligand_resmask).Selected())
    if not lig_sel:
        raise ValueError(f"ligand mask {ligand_resmask!r} selected no atoms")
    lig_heavy = [i for i in lig_sel if mol.atoms[i].element_name != "H"]
    if len(lig_heavy) < 3:
        raise ValueError("ligand has fewer than 3 heavy atoms; cannot build "
                         "a Boresch restraint")

    rec_sel = [i for i in pmd.amber.AmberMask(mol, receptor_mask).Selected()
               if i not in lig_sel]
    if len(rec_sel) < 3:
        raise ValueError(f"receptor mask {receptor_mask!r} gave <3 atoms")

    # ---- ligand COM (heavy atoms) ----
    com = np.mean([coords[i] for i in lig_heavy], axis=0)

    # ---- L1: heavy atom nearest COM, tie-broken by low RMSF ----
    lig_ranked = sorted(
        lig_heavy,
        key=lambda i: (np.linalg.norm(coords[i]-com), rmsf_of(i)))
    L1 = lig_ranked[0]

    # ---- candidate pools (kept small for a cheap exhaustive search) ----
    # L2/L3 from ligand heavy atoms; P1/P2/P3 from the nearest stable CA shell.
    shell = [i for i in rec_sel
             if REC_SHELL_MIN <= _dist(coords[L1], coords[i]) <= REC_SHELL_MAX]
    if len(shell) < 3:
        shell = sorted(rec_sel, key=lambda i: _dist(coords[L1], coords[i]))[:30]
    # rank CA candidates by (rmsf, distance) and keep a manageable top-K
    shell = sorted(shell, key=lambda i: (rmsf_of(i), _dist(coords[L1], coords[i])))[:12]
    lig_others = [i for i in lig_heavy if i != L1]

    def _angle_quality(a):
        """1.0 at 90 deg, 0.0 at the sanity edges, negative outside."""
        return 1.0 - abs(a - 90.0) / 60.0

    def _score(L2, L3, P1, P2, P3):
        # hard rejects
        if _dist(coords[L1], coords[L2]) < LIG_MIN_SEP: return None
        if _dist(coords[L1], coords[L3]) < LIG_MIN_SEP: return None
        if _dist(coords[L2], coords[L3]) < LIG_MIN_SEP: return None
        if _dist(coords[P1], coords[P2]) < REC_MIN_SEP: return None
        if _dist(coords[P1], coords[P3]) < REC_MIN_SEP: return None
        if _dist(coords[P2], coords[P3]) < REC_MIN_SEP: return None
        if _collinear(coords[L1], coords[L2], coords[L3]): return None
        if _collinear(coords[L1], coords[P1], coords[L2]): return None
        if _collinear(coords[P1], coords[P2], coords[P3]): return None
        alpha = _angle(coords[P1], coords[L1], coords[L2])  # P1-L1-L2
        theta = _angle(coords[P2], coords[P1], coords[L1])  # P2-P1-L1
        gamma_in = _angle(coords[L1], coords[L2], coords[L3])
        if not (SANITY_ANG_MIN <= alpha <= SANITY_ANG_MAX): return None
        if not (SANITY_ANG_MIN <= theta <= SANITY_ANG_MAX): return None
        # soft score: angle conditioning + low fluctuation + good spread
        rmsf_pen = -(rmsf_of(L2)+rmsf_of(L3)+rmsf_of(P1)+rmsf_of(P2)+rmsf_of(P3))
        return (_angle_quality(alpha) + _angle_quality(theta)
                + 0.5*_angle_quality(gamma_in)
                + 0.1*rmsf_pen)

    best = None; best_combo = None
    for P1 in shell:
        # only consider receptor pairs/triples within this P1's neighbourhood
        others = [i for i in shell if i != P1]
        for L2 in lig_others:
            for L3 in lig_others:
                if L3 == L2: continue
                for P2 in others:
                    for P3 in others:
                        if P3 == P2: continue
                        s = _score(L2, L3, P1, P2, P3)
                        if s is not None and (best is None or s > best):
                            best = s; best_combo = (L2, L3, P1, P2, P3)
    if best_combo is None:
        raise ValueError(
            "could not find a well-conditioned Boresch six-atom set "
            "automatically (ligand may be too small/floppy or the receptor "
            "shell too sparse). Provide manual masks via "
            "lig_restraint_atoms / rec_restraint_atoms.")
    L2, L3, P1, P2, P3 = best_combo

    lig_atoms = [L1, L2, L3]
    rec_atoms = [P1, P2, P3]
    lig_masks = [amber_idx(mol.atoms[i]) for i in lig_atoms]
    rec_masks = [amber_idx(mol.atoms[i]) for i in rec_atoms]

    # ---- geometry + sanity report ----
    checks = _validate(coords, lig_atoms, rec_atoms, rmsf_of)
    if report_path:
        _write_report(report_path, mol, lig_atoms, rec_atoms,
                      lig_masks, rec_masks, coords, rmsf_of, checks)
    if not checks["all_pass"]:
        raise ValueError(
            "auto-selected Boresch atoms failed geometric sanity "
            f"(see report): {checks['failures']}. Provide manual masks.")
    return lig_masks, rec_masks

def _validate(coords, lig, rec, rmsf_of):
    L1, L2, L3 = lig
    P1, P2, P3 = rec
    r = _dist(coords[L1], coords[P1])
    alpha = _angle(coords[P1], coords[L1], coords[L2])    # P1-L1-L2
    theta = _angle(coords[P2], coords[P1], coords[L1])    # P2-P1-L1
    failures = []
    if not (REC_SHELL_MIN - 1 <= r <= REC_SHELL_MAX + 2):
        failures.append(f"r={r:.2f} out of range")
    for nm, a in (("alpha", alpha), ("theta", theta)):
        if not (SANITY_ANG_MIN <= a <= SANITY_ANG_MAX):
            failures.append(f"{nm}={a:.1f} outside [{SANITY_ANG_MIN},{SANITY_ANG_MAX}]")
    if _collinear(coords[P1], coords[L1], coords[L2]):
        failures.append("P1-L1-L2 collinear")
    if _collinear(coords[P2], coords[P1], coords[L1]):
        failures.append("P2-P1-L1 collinear")
    return dict(r=r, alpha=alpha, theta=theta,
                failures=failures, all_pass=(len(failures) == 0))


def _load_parm_with_coords(prmtop, inpcrd):
    """Load prmtop + coords robustly. v2.5.31h: Amber-22 parmed calls
    np.array(box, copy=False) which raises under NumPy 2.0 -> the post-eq Boresch
    gate crashed (rc=1, production refused). Fall back to loading topology and
    coordinates separately; else raise a clear 'pip install numpy<2' message."""
    import parmed as pmd
    try:
        return pmd.load_file(str(prmtop), xyz=str(inpcrd))
    except (ValueError, TypeError) as e:
        msg = str(e)
        is_numpy2 = ("Unable to avoid copy" in msg or "copy=False" in msg
                     or "copy keyword" in msg)
        try:
            parm = pmd.load_file(str(prmtop))
            rst = pmd.load_file(str(inpcrd))
            parm.coordinates = rst.coordinates
            try:
                parm.box = rst.box
            except Exception:
                pass
            return parm
        except Exception:
            if is_numpy2:
                raise RuntimeError(
                    "parmed failed because NumPy 2.x is on the path but Amber-22 "
                    "parmed needs NumPy 1.x. Fix on the node: pip install 'numpy<2' "
                    "then pip install -e . --no-deps. Original error: " + msg) from e
            raise


def validate_masks(prmtop, inpcrd, lig_masks, rec_masks, lig_idx=0):
    """Validate an EXISTING six-atom set against a structure (e.g. the
    equilibrated rst7). Returns the same dict as _validate() plus the resolved
    0-based indices. Used by the post-equilibration fail-fast gate.
    """
    import parmed as pmd
    mol = _load_parm_with_coords(prmtop, inpcrd)
    coords = mol.coordinates

    def _resolve(mask):
        sel = list(pmd.amber.AmberMask(mol, mask).Selected())
        if not sel:
            raise ValueError(f"mask {mask!r} selected no atoms")
        return sel[lig_idx if lig_idx < len(sel) else 0]

    lig = [_resolve(m) for m in lig_masks]
    rec = [_resolve(m) for m in rec_masks]
    checks = _validate(coords, lig, rec, lambda i: 0.0)
    checks["lig_idx0"] = lig
    checks["rec_idx0"] = rec
    return checks

def _write_report(path, mol, lig, rec, lig_masks, rec_masks, coords,
                  rmsf_of, checks):
    L = ["# Auto-selected Boresch six-atom restraint (amber_md v2.5.17)",
         "# RMSF-aware selection; Calpha receptor anchors.",
         "",
         "role  mask              atom            rmsf(A)"]
    roles = ["L1", "L2", "L3", "P1", "P2", "P3"]
    for role, idx, mask in zip(roles, lig + rec, lig_masks + rec_masks):
        at = mol.atoms[idx]
        L.append(f"{role:4s}  {mask:16s}  {at.residue.name}{at.residue.idx+1}/"
                 f"{at.name:4s}  {rmsf_of(idx):7.3f}")
    L += ["",
          f"r (L1-P1)     = {checks['r']:.3f} A",
          f"alpha(P1-L1-L2)= {checks['alpha']:.2f} deg",
          f"theta(P2-P1-L1)= {checks['theta']:.2f} deg",
          "",
          "sanity: " + ("PASS" if checks["all_pass"]
                        else "FAIL -> " + "; ".join(checks["failures"])),
          "",
          "lig_restraint_atoms: " + str(lig_masks),
          "rec_restraint_atoms: " + str(rec_masks)]
    Path(path).write_text("\n".join(L) + "\n")


def select_or_manual(prmtop, inpcrd, lig_restraint_atoms=None,
                     rec_restraint_atoms=None, **kw):
    """Convenience wrapper used by the driver.

    If BOTH lig_restraint_atoms and rec_restraint_atoms are 3-element lists,
    they are returned unchanged (manual override, matching FEP-SPell). Else
    autoselect_boresch_atoms() runs.
    """
    manual_lig = isinstance(lig_restraint_atoms, (list, tuple)) and len(lig_restraint_atoms) == 3
    manual_rec = isinstance(rec_restraint_atoms, (list, tuple)) and len(rec_restraint_atoms) == 3
    if manual_lig and manual_rec:
        return list(lig_restraint_atoms), list(rec_restraint_atoms)
    return autoselect_boresch_atoms(prmtop, inpcrd, **kw)


if __name__ == "__main__":
    # geometry self-tests (no parmed needed): exercise the vector helpers
    A = (0.0, 0.0, 0.0); B = (1.0, 0.0, 0.0); C = (1.0, 1.0, 0.0)
    assert abs(_dist(A, B) - 1.0) < 1e-9
    assert abs(_angle(A, B, C) - 90.0) < 1e-6, _angle(A, B, C)
    assert _collinear((0, 0, 0), (1, 0, 0), (2, 0, 0))
    assert not _collinear((0, 0, 0), (1, 0, 0), (1, 1, 0))
    # validate() on a clean geometry passes; a collinear one fails
    coords = {0:(0,0,0),1:(0.0,1.5,0),2:(0.0,1.5,1.5),  # L1,L2,L3 (L2 off the L1-P1 axis)
              3:(7.0,0,0),4:(7.0,5.0,0),5:(11.0,5.0,3.0)}  # P1,P2,P3
    ok = _validate(coords, [0,1,2], [3,4,5], lambda i:0.0)
    print("clean geometry checks:", ok["all_pass"], ok["failures"])
    assert ok["all_pass"], ok
    bad = _validate({0:(0,0,0),1:(1.5,0,0),2:(3,0,0),3:(7,0,0),4:(14,0,0),5:(21,0,0)},
                    [0,1,2],[3,4,5], lambda i:0.0)
    print("collinear geometry checks:", bad["all_pass"], bad["failures"])
    assert not bad["all_pass"]
    print("\nself-test OK")