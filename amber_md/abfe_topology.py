"""ABFE topology builders.

v2.4.21:
  * NEW: GPU small-box guard in build_solvent_only_topology(). After
    computing the solvation buffer, estimate the truncated-octahedron
    perpendicular width and auto-raise the buffer if it would fall
    below 3*(cutoff+skinnb) (+10%% NPT margin) -- the pmemd.cuda
    minimum. Prevents the rc=255 small-box crash at decoupled vdw/
    decharge endpoints that a 12 A buffer (=> ~42.5 A box) caused.
    Per Amber dev guidance, -AllowSmallBox is unsafe; a bigger box is
    the correct fix.

v2.4.13 hotfix:
  * build_solvent_only_topology(): fix tleap exit 21 when ligand is neutral
    and ligand-only system needs solvation.
    - Skip the redundant "addIonsRand system Na+ 0 Cl- 0" line when the
      ligand has no net charge (tleap >=22 errors on it).
    - Use addIons2 (grid-based, more robust on small boxes) instead of
      addIonsRand when only adding salt to a single-residue system.
    - Use the ligand UNIT directly for solvation (skip `combine`) so
      tleap gets a proper bounding box. `combine` of a single unit drops
      coordinates in some tleap versions.
    - Add `check system` before save so any remaining issue is loud.
    - Add `loadOff atomic_ions` if not auto-loaded, for older Amber that
      doesn't pull ion types from leaprc.water.tip3p alone.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .config import SystemConfig
from .logger import get_logger
from .utils import ensure_dir, run, which_or_die, CommandError
from .prep import (PDBCleaner, LigandParametrizer,
                   _count_sdf_records, _count_mol2_records)
from .protonation import apply_protonation
from .builder import SystemBuilder, _WATER_MAP

log = get_logger()


def _extract_sdf_record(src: Path, index: int, dest: Path) -> Path:
    seen = 0
    buf: list[str] = []
    with open(src, "r", errors="replace") as f:
        for line in f:
            buf.append(line)
            if line.strip() == "$$$$":
                if seen == index:
                    dest.write_text("".join(buf))
                    return dest
                seen += 1
                buf = []
    if seen == index and buf:
        if not buf[-1].rstrip().endswith("$$$$"):
            buf.append("$$$$\n")
        dest.write_text("".join(buf))
        return dest
    raise IndexError(
        f"SDF/MOL {src} has at most {max(seen, 0)} record(s); "
        f"asked for index {index}.")


def _extract_mol2_record(src: Path, index: int, dest: Path) -> Path:
    records: list[list[str]] = []
    current: list[str] = []
    with open(src, "r", errors="replace") as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                if current:
                    records.append(current)
                current = [line]
            else:
                if current:
                    current.append(line)
        if current:
            records.append(current)
    if index >= len(records):
        raise IndexError(
            f"MOL2 {src} has {len(records)} record(s); asked for index {index}.")
    dest.write_text("".join(records[index]))
    return dest


def extract_ligand_record(ligand_file: Path, index: int, dest_dir: Path) -> Path:
    ligand_file = Path(ligand_file).resolve()
    ext = ligand_file.suffix.lower()
    dest_dir = ensure_dir(dest_dir)
    if ext in (".sdf", ".mol"):
        n = _count_sdf_records(ligand_file)
    elif ext == ".mol2":
        n = _count_mol2_records(ligand_file)
    else:
        n = 1
    if n <= 1 and index == 0:
        log.info("Ligand %s holds 1 record; using as-is.", ligand_file.name)
        return ligand_file
    if index >= n and n > 0:
        raise IndexError(
            f"{ligand_file} has {n} record(s); asked for index {index}.")
    out = dest_dir / f"ligand_{index}{ext}"
    log.info("Extracting record %d from %s -> %s", index, ligand_file, out)
    if ext in (".sdf", ".mol"):
        return _extract_sdf_record(ligand_file, index, out)
    if ext == ".mol2":
        return _extract_mol2_record(ligand_file, index, out)
    raise ValueError(
        f"Cannot extract record {index} from {ext} file; "
        f"only .sdf/.mol/.mol2 support multi-record splitting.")


def build_abfe_topology(
    protein_pdb: Path,
    ligand_file: Path,
    ligand_index: int,
    ligand_resname: str,
    work_dir: Path,
    sys_cfg: Optional[SystemConfig] = None,
    auto_protonation: bool = True,
    protonation_overrides: Optional[dict] = None,
) -> tuple[Path, Path, Path]:
    if sys_cfg is None:
        sys_cfg = SystemConfig()

    work_dir = Path(work_dir).expanduser().resolve()
    ensure_dir(work_dir)

    one_lig = extract_ligand_record(
        Path(ligand_file).resolve(), ligand_index, work_dir / "ligand_input")

    cleaner = PDBCleaner(work_dir / "prep")
    log.info("ABFE topology: cleaning protein %s", protein_pdb)
    prot = cleaner.clean_protein_only(Path(protein_pdb).resolve())
    if auto_protonation:
        prot_fixed = prot.parent / (prot.stem + "_proton.pdb")
        log.info("ABFE topology: applying protonation -> %s", prot_fixed)
        apply_protonation(prot, prot_fixed,
                          manual_overrides=protonation_overrides)
        prot = prot_fixed

    lp = LigandParametrizer(work_dir, sys_cfg)
    log.info("ABFE topology: parametrizing ligand (resname=%s)", ligand_resname)
    mol2, frcmod = lp.parametrize(one_lig, ligand_resname)

    builder = SystemBuilder(work_dir, sys_cfg)
    log.info("ABFE topology: building solvated complex (tleap)")
    prmtop, inpcrd, complex_pdb = builder.build(
        prot, mol2, frcmod, ligand_resname)

    import json
    (work_dir / "abfe_topology_inputs.json").write_text(json.dumps({
        "protein_pdb":   str(protein_pdb),
        "ligand_file":   str(ligand_file),
        "ligand_index":  ligand_index,
        "ligand_resname": ligand_resname,
        "prmtop":        str(prmtop),
        "inpcrd":        str(inpcrd),
        "complex_pdb":   str(complex_pdb),
        "ligand_mol2":   str(mol2),
        "ligand_frcmod": str(frcmod),
    }, indent=2))

    log.info("ABFE topology done: prmtop=%s  inpcrd=%s  complex_pdb=%s",
             prmtop, inpcrd, complex_pdb)
    return Path(prmtop), Path(inpcrd), Path(complex_pdb)


# =====================================================================
# v2.4.13 — solvent-leg topology builder (HOTFIX for tleap exit 21)
# =====================================================================
def build_solvent_only_topology(
    ligand_file: Path,
    ligand_index: int,
    ligand_resname: str,
    work_dir: Path,
    sys_cfg: Optional[SystemConfig] = None,
    reuse_mol2: Optional[Path] = None,
    reuse_frcmod: Optional[Path] = None,
) -> tuple[Path, Path]:
    """Build a ligand-in-water box for the ABFE solvent leg.

    v2.4.13 fixes for tleap exit 21:
      * Don't emit "addIonsRand system Na+ 0 Cl- 0" when ligand is neutral
        (tleap 22+ rejects it).
      * Use the ligand UNIT directly (skip `combine`) so the bounding box
        is preserved for solvateOct/solvateBox.
      * Use addIons2 for salt addition (grid-based; robust on small boxes).
    """
    if sys_cfg is None:
        sys_cfg = SystemConfig()

    work_dir = Path(work_dir).expanduser().resolve()
    ensure_dir(work_dir)

    if reuse_mol2 is not None and reuse_frcmod is not None:
        mol2 = Path(reuse_mol2).resolve()
        frcmod = Path(reuse_frcmod).resolve()
        if not mol2.exists():
            raise FileNotFoundError(f"reuse_mol2 not found: {mol2}")
        if not frcmod.exists():
            raise FileNotFoundError(f"reuse_frcmod not found: {frcmod}")
        log.info("Solvent leg: reusing ligand params from complex leg")
        log.info("  mol2   = %s", mol2)
        log.info("  frcmod = %s", frcmod)
    else:
        log.warning("Solvent leg: re-parametrizing ligand from scratch.")
        one_lig = extract_ligand_record(
            Path(ligand_file).resolve(),
            ligand_index, work_dir / "ligand_input")
        lp = LigandParametrizer(work_dir, sys_cfg)
        mol2, frcmod = lp.parametrize(one_lig, ligand_resname)

    which_or_die("tleap")
    build_dir = ensure_dir(work_dir / "build_solvent")
    prmtop = build_dir / "solvent.prmtop"
    inpcrd = build_dir / "solvent.inpcrd"
    outpdb = build_dir / "solvent.pdb"
    leap_in = build_dir / "tleap_solvent.in"

    water_leaprc, box = _WATER_MAP[sys_cfg.water_model.lower()]
    solv_cmd = ("solvateOct" if sys_cfg.box_shape == "octahedral"
                else "solvateBox")
    buf = sys_cfg.box_buffer_A

    # --- v2.4.21: GPU small-box guard ----------------------------------
    # pmemd.cuda requires >= 3 hash cells per dimension; a hash cell is
    # ~(cutoff + skinnb). For a truncated octahedron the perpendicular
    # width is ~0.77 * the tleap "iso" edge, and NPT shrinks the box
    # further. A 12 A buffer on a small ligand gave a ~42.5 A box that
    # dropped below the limit at decoupled endpoints (rc=255 crash).
    # Auto-raise the buffer so min_perp_width comfortably clears the
    # threshold. See Amber dev (D. Case): -AllowSmallBox is NOT safe;
    # the correct fix is a bigger box.
    try:
        _cut  = float(getattr(sys_cfg, "cutoff_A", 10.0) or 10.0)
        _skin = float(getattr(sys_cfg, "skinnb_A", 2.0) or 2.0)
        _is_oct = (sys_cfg.box_shape == "octahedral")
        # need >= 3 cells per perpendicular dim, with NPT-shrink margin
        _min_perp_needed = 3.0 * (_cut + _skin) * 1.10   # 10% NPT margin
        # tleap "iso" edge ≈ solute_max_dim + 2*buf; perp width ≈ 0.77*edge (oct)
        # Conservatively assume a near-point solute (worst case): edge ≈ 2*buf.
        _perp_factor = 0.77 if _is_oct else 1.0
        _est_perp = _perp_factor * (2.0 * buf)
        if _est_perp < _min_perp_needed:
            _new_buf = (_min_perp_needed / _perp_factor) / 2.0
            _new_buf = max(_new_buf, buf)  # never shrink
            log.warning(
                "GPU small-box guard: buffer %.1f A -> est. perp width %.1f A "
                "< required %.1f A (cut=%.1f, skin=%.1f, oct=%s). "
                "Auto-raising buffer to %.1f A to avoid pmemd.cuda small-box "
                "crash at decoupled endpoints.",
                buf, _est_perp, _min_perp_needed, _cut, _skin, _is_oct, _new_buf)
            buf = round(_new_buf + 0.5, 1)  # round up to 0.1 A
    except Exception as _e:
        log.warning("GPU small-box guard skipped (%s); using buffer=%.1f A",
                    _e, buf)
    # -------------------------------------------------------------------

    # v2.4.13: ligand charge defaults to 0; only add neutralizing ions
    # if the ligand actually carries a net charge. Always honor user-
    # requested salt concentration.
    lig_charge = int(getattr(sys_cfg, "ligand_charge", 0) or 0)

    ion_lines = []
    if sys_cfg.neutralize and lig_charge != 0:
        # addIons2 picks the counter-ion based on system net charge sign
        if lig_charge > 0:
            ion_lines.append(f"addIons2 system Cl- {lig_charge}")
        else:
            ion_lines.append(f"addIons2 system Na+ {-lig_charge}")
    if sys_cfg.salt_conc_M > 0:
        n_salt = int(round(sys_cfg.salt_conc_M * 100))
        if n_salt > 0:
            # addIons2 (grid) is more robust than addIonsRand on small boxes
            ion_lines.append(f"addIons2 system Na+ {n_salt} Cl- {n_salt}")
    ion_block = "\n".join(ion_lines) if ion_lines else "# no ions needed"

    # v2.4.13: solvate the ligand UNIT directly, no `combine` wrapper.
    leap_in.write_text(f"""# Solvent-leg topology: ligand-in-water only.
# Generated by amber_md.abfe_topology.build_solvent_only_topology (v2.4.13)
verbosity 1
source leaprc.{sys_cfg.ligand_ff}
source {water_leaprc}
loadAmberParams {frcmod}
{ligand_resname} = loadMol2 {mol2}
# Rename the unit handle to 'system' so the rest of the script is
# identical to the complex pathway.
system = {ligand_resname}
{solv_cmd} system {box} {buf}
{ion_block}
check system
saveAmberParm system {prmtop} {inpcrd}
savePDB system {outpdb}
quit
""")
    try:
        run(["tleap", "-f", str(leap_in)], cwd=build_dir)
    except CommandError as e:
        # tleap writes its real diagnostics to leap.log, not stderr.
        log_file = build_dir / "leap.log"
        if log_file.exists():
            tail = "\n".join(log_file.read_text().splitlines()[-40:])
            raise CommandError(
                f"{e}\n--- tail of {log_file} ---\n{tail}\n--- end leap.log ---"
            ) from e
        raise

    if not prmtop.exists():
        log_file = build_dir / "leap.log"
        tail = (log_file.read_text().splitlines()[-40:]
                if log_file.exists() else ["(no leap.log)"])
        raise RuntimeError(
            "tleap finished without producing solvent.prmtop. Tail of leap.log:\n"
            + "\n".join(tail))

    import json
    (build_dir / "solvent_topology_inputs.json").write_text(json.dumps({
        "ligand_file":    str(ligand_file),
        "ligand_index":   ligand_index,
        "ligand_resname": ligand_resname,
        "ligand_mol2":    str(mol2),
        "ligand_frcmod":  str(frcmod),
        "prmtop":         str(prmtop),
        "inpcrd":         str(inpcrd),
        "ligand_charge":  lig_charge,
        "reused_params":  reuse_mol2 is not None,
    }, indent=2))

    log.info("Solvent leg topology done: prmtop=%s  inpcrd=%s",
             prmtop, inpcrd)
    return Path(prmtop), Path(inpcrd)
