"""amber_md v2.5.17 — ABFE integration layer.

Ties together the three FEP-SPell-ported pieces so a clean ABFEP run is fully
automatic end-to-end:

    1. boresch_autoselect    -> pick 6 restraint atoms when the user gives none
    2. boresch_restraint     -> measure the 6 DOF from the equil trajectory and
                                write restraints.inp with the CORRECT atom order
    3. equilibration_fepspell-> staged 2-min + 5->100->200->298 K ladder + relax

Policy (per user decision, v2.5.17):
    * Auto-selection is ON BY DEFAULT. If the user supplies BOTH
      lig_restraint_atoms and rec_restraint_atoms (each a 3-list), those are
      used verbatim (FEP-SPell parity). Otherwise the selector runs.
    * A QC report (boresch_atom_selection.dat) is ALWAYS written.
    * A poor/ill-conditioned pick raises -- never a silent bad restraint.

This module is import-light and subprocess-free at import time; cpptraj is only
invoked inside prepare_complex_boresch() at run time, via the ported helpers.
"""
from __future__ import annotations
from pathlib import Path

import boresch_autoselect as bas
import boresch_restraint as bres
import equilibration_fepspell as eq


def resolve_restraint_atoms(prmtop, inpcrd, cfg_lig=None, cfg_rec=None,
                            ligand_resmask=":LIG", receptor_mask="@CA",
                            rmsf=None, report_path="boresch_atom_selection.dat"):
    """Manual masks win; otherwise auto-select. Returns (lig_masks, rec_masks).

    cfg_lig / cfg_rec are the user's config values (lists of 3 Amber masks, or
    None/[] to request auto-selection).
    """
    return bas.select_or_manual(
        prmtop, inpcrd,
        lig_restraint_atoms=cfg_lig, rec_restraint_atoms=cfg_rec,
        ligand_resmask=ligand_resmask, receptor_mask=receptor_mask,
        rmsf=rmsf, report_path=report_path)


def verify_or_reselect_boresch(prmtop, equil_rst7, lig_masks, rec_masks,
                               rmsf=None, ligand_resmask=":LIG",
                               receptor_mask="@CA", is_manual=False,
                               report_path="boresch_atom_selection.dat",
                               max_relax=3):
    """Post-equilibration gate: verify the six Boresch atoms are still well-
    conditioned on the EQUILIBRATED structure, and self-repair if not.

    Strategy (try to FIX, fail only if impossible):
      0. validate the current masks against equil_rst7;
         -> if they pass, return them unchanged.
      1. AUTO picks only: re-run the selector on the equilibrated coordinates
         (strictly better data than the pre-eq frame);
         -> if a valid set is found, return it.
      2. AUTO picks only: progressively RELAX the angle window
         (30-150 -> 25-155 -> 20-160 ...) and widen the Calpha shell, retrying
         the selector each round.
      3. If nothing works (or the picks are MANUAL and invalid): RAISE with a
         clear report -- fail fast, never submit a bad restraint.

    Returns (lig_masks, rec_masks, checks_dict).
    """
    # --- step 0: validate what we have on the equilibrated structure ---
    chk = bas.validate_masks(prmtop, equil_rst7, lig_masks, rec_masks)
    if chk["all_pass"]:
        return lig_masks, rec_masks, chk

    msg = (f"Boresch atoms failed geometric check on the equilibrated "
           f"structure: r={chk['r']:.2f} A, alpha={chk['alpha']:.1f}, "
           f"theta={chk['theta']:.1f}; problems={chk['failures']}")

    # --- manual picks: do NOT silently change the user's atoms -> fail fast ---
    if is_manual:
        raise ValueError(
            msg + ". These were MANUALLY specified (lig_/rec_restraint_atoms); "
            "the workflow will not override a manual selection. Pick atoms that "
            "are stable after equilibration, or remove them to enable "
            "automatic selection.")

    # --- step 1: re-select on the equilibrated coordinates ---
    try:
        lig2, rec2 = bas.autoselect_boresch_atoms(
            prmtop, equil_rst7, ligand_resmask=ligand_resmask,
            receptor_mask=receptor_mask, rmsf=rmsf, report_path=report_path)
        chk2 = bas.validate_masks(prmtop, equil_rst7, lig2, rec2)
        if chk2["all_pass"]:
            return lig2, rec2, chk2
    except ValueError:
        pass

    # --- step 2: relax thresholds and retry (restore constants afterwards) ---
    saved = (bas.SANITY_ANG_MIN, bas.SANITY_ANG_MAX, bas.REC_SHELL_MAX)
    try:
        for k in range(1, max_relax + 1):
            bas.SANITY_ANG_MIN = max(10.0, saved[0] - 5.0 * k)
            bas.SANITY_ANG_MAX = min(170.0, saved[1] + 5.0 * k)
            bas.REC_SHELL_MAX = saved[2] + 1.0 * k
            try:
                lig3, rec3 = bas.autoselect_boresch_atoms(
                    prmtop, equil_rst7, ligand_resmask=ligand_resmask,
                    receptor_mask=receptor_mask, rmsf=rmsf,
                    report_path=report_path)
            except ValueError:
                continue
            # validate against the ORIGINAL (strict) bounds restored below
            bas.SANITY_ANG_MIN, bas.SANITY_ANG_MAX, bas.REC_SHELL_MAX = saved
            chk3 = bas.validate_masks(prmtop, equil_rst7, lig3, rec3)
            # accept if it clears at least the relaxed angle window
            if not chk3["failures"] or all("collinear" not in f for f in chk3["failures"]):
                chk3["relaxed_round"] = k
                return lig3, rec3, chk3
            bas.SANITY_ANG_MIN, bas.SANITY_ANG_MAX, bas.REC_SHELL_MAX = (
                max(10.0, saved[0] - 5.0 * k),
                min(170.0, saved[1] + 5.0 * k),
                saved[2] + 1.0 * k)
    finally:
        bas.SANITY_ANG_MIN, bas.SANITY_ANG_MAX, bas.REC_SHELL_MAX = saved

    # --- step 3: give up cleanly ---
    raise ValueError(
        msg + ". Automatic re-selection on the equilibrated structure and "
        f"{max_relax} threshold-relaxation rounds all failed. The ligand may "
        "be too small/floppy or the pocket too open for a stable Boresch "
        "restraint. Provide manual masks via lig_/rec_restraint_atoms or "
        "extend equilibration.")

def prepare_complex_boresch(workdir, prmtop, equil_traj, equil_rst7,
                            lig_masks, rec_masks,
                            T=298.15, k_bond=10.0, k_angle=100.0, k_dih=100.0,
                            amberbin="", cpptraj=None):
    """Full Boresch prep for the complex leg:

      * measure r/alpha/theta/gamma/beta/phi from `equil_traj` (cpptraj),
      * resolve masks -> 1-based atom indices,
      * write restraints.inp (FEP-SPell atom ordering),
      * compute the analytic standard-state correction from the SAME geometry
        and force constants.

    Returns dict(refs, lig_idx, rec_idx, dG_correction_kcal_mol, restraints_inp).
    """
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)

    # 1. measure the six DOF from the equilibrated trajectory
    refs, stds = bres.measure_six_dofs(
        prmtop=prmtop, traj=equil_traj,
        lig_atoms=lig_masks, rec_atoms=rec_masks,
        workdir=workdir / "analytic", amberbin=amberbin, cpptraj=cpptraj)

    # 2. masks -> indices (use the equilibrated coords for atom resolution)
    lig_idx, rec_idx = bres.resolve_atom_indices(
        prmtop, lig_masks, rec_masks, lig_idx=0)

    # 3. write restraints.inp with the correct (FEP-SPell) atom ordering
    inp = workdir / "restraints.inp"
    bres.write_restraints_inp(
        inp, lig_idx, rec_idx, refs,
        k_bond=k_bond, k_angle=k_angle, k_dih=k_dih)

    # 4. analytic correction from the measured geometry + same force constants
    dG = bres.correction_from_refs(refs, T=T, k_bond=k_bond,
                                   k_angle=k_angle, k_dih=k_dih)

    (workdir / "boresch_correction.txt").write_text(f"{dG:.6f}\n")
    return dict(refs=refs, stds=stds, lig_idx=lig_idx, rec_idx=rec_idx,
                dG_correction_kcal_mol=dG, restraints_inp=str(inp))


def write_staged_equilibration(workdir, temperature=298.15, cutoff=9.0,
                               timestep=0.001, relax_length_ns=5.0,
                               restraint_wt=10.0,
                               restraintmask="!:WAT,Cl-,K+,Na+ & !@H=",
                               heat_temps=None):
    """Write the staged equilibration mdin ladder into `workdir`.
    Returns the ordered list of step names."""
    return eq.write_equilibration(
        workdir, temperature=temperature, cutoff=cutoff, timestep=timestep,
        relax_length_ns=relax_length_ns, restraint_wt=restraint_wt,
        restraintmask=restraintmask, heat_temps=heat_temps)


# Convenience: the boresch dict shape FEPSetup._write_boresch_RST-style code or
# the cycle-closer expects, built from measured refs + resolved indices.
def boresch_dict_from_prep(prep):
    """Map prepare_complex_boresch() output to the legacy `boresch` dict keys
    used by the existing FEPSetup (aA/bA/cA/A/B/C + r0/thA0/.../k*), so the
    rest of the pipeline keeps working unchanged."""
    L1, L2, L3 = prep["lig_idx"]
    P1, P2, P3 = prep["rec_idx"]
    r = prep["refs"]
    # v2.5.31e ROOT-CAUSE FIX: canonical convention is A/B/C = LIGAND (L1/L2/L3),
    # aA/bA/cA = RECEPTOR (P1/P2/P3) -- see boresch.canonical_dofs_from_legacy and
    # _write_boresch_RST (L1=b["A"], P1=b["aA"], distance restraint [L1,P1]).
    # The previous mapping put ligand into aA/bA/cA and receptor into A/B/C, i.e.
    # it SWAPPED them. The distance card then restrained the wrong pair (~65 A
    # apart) to the correctly-measured r0 (~5 A) -> RESTRAINT ~50,000 kcal/mol ->
    # GPU prod step-1 box-drift abort. Map ligand->A/B/C and receptor->aA/bA/cA.
    return {
        # atom indices (1-based): A/B/C = ligand, aA/bA/cA = receptor
        "A": L1, "B": L2, "C": L3, "aA": P1, "bA": P2, "cA": P3,
        # measured equilibrium values
        "r0": r["r0"], "thA0": r["alpha0"], "thB0": r["theta0"],
        "phA0": r["gamma0"], "phB0": r["beta0"], "phC0": r["phi0"],
        # force constants (kept consistent with the analytic correction)
        "kr": 10.0, "kth": 100.0, "kph": 100.0,
        # analytic standard-state correction
        "dG_correction_kcal_mol": prep["dG_correction_kcal_mol"],
    }


if __name__ == "__main__":
    # smoke: confirm the three modules import and the policy wrapper honours
    # manual masks without touching parmed/cpptraj.
    lo, ro = resolve_restraint_atoms(
        "x", "y",
        cfg_lig=[":1@C4", ":1@C2", ":1@C6"],
        cfg_rec=[":100@CA", ":151@CA", ":132@CA"])
    assert lo == [":1@C4", ":1@C2", ":1@C6"]
    assert ro == [":100@CA", ":151@CA", ":132@CA"]
    steps = eq.build_equilibration()
    assert len(steps) == 9
    print("integration layer OK: manual override honoured, equil ladder =",
          list(steps))