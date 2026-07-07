# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.
"""Per-edge dual-topology builder for Amber RBFE (Option B, final58).

For a relative free-energy edge A->B we build TWO solvated systems in which
*both* ligands coexist so pmemd TI (timask1=:L1, timask2=:L2) can morph A->B:
  * complex : protein + ligand_A(:L1) + ligand_B(:L2) in water
  * solvent : ligand_A(:L1) + ligand_B(:L2) in water (no protein)
Reuses the tested prep stack (PDBCleaner, LigandParametrizer, SystemConfig,
_WATER_MAP) and only adds the tleap glue that combines two ligand units.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import json

from .config import SystemConfig
from .prep import PDBCleaner, LigandParametrizer
from .protonation import apply_protonation
from .abfe_topology import extract_ligand_record
from .builder import _WATER_MAP
from .utils import run, which_or_die, ensure_dir
from .logger import get_logger

log = get_logger()

RESNAME_A = "L1"
RESNAME_B = "L2"


def _ion_block(cfg: SystemConfig) -> str:
    lines = []
    n_salt = int(round(cfg.salt_conc_M * 100)) if cfg.salt_conc_M > 0 else 0
    if cfg.neutralize:
        lines.append("addIons2 system Na+ 0")
        lines.append("addIons2 system Cl- 0")
    if n_salt > 0:
        lines.append(f"addIonsRand system Na+ {n_salt} Cl- {n_salt}")
    return "\n".join(lines) if lines else "# no ions needed"


def _solvate_cmd(cfg: SystemConfig, buf: float) -> str:
    _, box = _WATER_MAP[cfg.water_model.lower()]
    kw = "solvateOct" if cfg.box_shape == "octahedral" else "solvateBox"
    return f"{kw} system {box} {buf}"


def _parametrize_two(lig_a_file: Path, lig_b_file: Path,
                     work_dir: Path, cfg: SystemConfig):
    lp_a = LigandParametrizer(ensure_dir(work_dir / "ligand_A"), cfg)
    lp_b = LigandParametrizer(ensure_dir(work_dir / "ligand_B"), cfg)
    rec_a = extract_ligand_record(Path(lig_a_file).resolve(), 0,
                                  work_dir / "ligand_A_input")
    rec_b = extract_ligand_record(Path(lig_b_file).resolve(), 0,
                                  work_dir / "ligand_B_input")
    log.info("RBFE edge: parametrizing ligand A (resname=%s)", RESNAME_A)
    mol2_a, frc_a = lp_a.parametrize(rec_a, RESNAME_A)
    log.info("RBFE edge: parametrizing ligand B (resname=%s)", RESNAME_B)
    mol2_b, frc_b = lp_b.parametrize(rec_b, RESNAME_B)
    return Path(mol2_a), Path(frc_a), Path(mol2_b), Path(frc_b)


def _tleap_build(leap_in: Path, build_dir: Path, prmtop: Path):
    run(["tleap", "-f", str(leap_in)], cwd=build_dir)
    if not prmtop.exists():
        leaplog = build_dir / "leap.log"
        tail = ("\n".join(leaplog.read_text().splitlines()[-40:])
                if leaplog.exists() else "(no leap.log)")
        raise RuntimeError(
            f"tleap failed to produce {prmtop}.\n--- tail of leap.log ---\n"
            f"{tail}\n--- end leap.log ---")


def build_rbfe_edge_topology(
    protein_pdb: Path,
    ligand_a_file: Path,
    ligand_b_file: Path,
    work_dir: Path,
    sys_cfg: Optional[SystemConfig] = None,
    auto_protonation: bool = True,
    protonation_overrides: Optional[dict] = None,
) -> dict:
    if sys_cfg is None:
        sys_cfg = SystemConfig()
    work_dir = Path(work_dir).expanduser().resolve()
    ensure_dir(work_dir)
    which_or_die("tleap")

    cleaner = PDBCleaner(work_dir / "prep")
    log.info("RBFE edge: cleaning protein %s", protein_pdb)
    prot = cleaner.clean_protein_only(Path(protein_pdb).resolve())
    if auto_protonation:
        prot_fixed = prot.parent / (prot.stem + "_proton.pdb")
        log.info("RBFE edge: applying protonation -> %s", prot_fixed)
        apply_protonation(prot, prot_fixed, manual_overrides=protonation_overrides)
        prot = prot_fixed

    mol2_a, frc_a, mol2_b, frc_b = _parametrize_two(
        Path(ligand_a_file), Path(ligand_b_file), work_dir, sys_cfg)

    water_leaprc, _ = _WATER_MAP[sys_cfg.water_model.lower()]
    buf = sys_cfg.box_buffer_A

    cbuild = ensure_dir(work_dir / "build")
    c_prm = cbuild / "complex.prmtop"; c_crd = cbuild / "complex.inpcrd"
    c_pdb = cbuild / "complex_solv.pdb"; c_leap = cbuild / "tleap.in"
    c_leap.write_text(f"""# RBFE dual-topology COMPLEX leg (A->B).
source leaprc.protein.{sys_cfg.protein_ff}
source leaprc.{sys_cfg.ligand_ff}
source {water_leaprc}
loadAmberParams {frc_a.resolve()}
loadAmberParams {frc_b.resolve()}
{RESNAME_A} = loadMol2 {mol2_a.resolve()}
{RESNAME_B} = loadMol2 {mol2_b.resolve()}
protein = loadPDB {prot.resolve()}
system = combine {{ protein {RESNAME_A} {RESNAME_B} }}
{_solvate_cmd(sys_cfg, buf)}
{_ion_block(sys_cfg)}
saveAmberParm system {c_prm} {c_crd}
savePDB system {c_pdb}
quit
""")
    log.info("RBFE edge: building solvated dual-topology complex (tleap)")
    _tleap_build(c_leap, cbuild, c_prm)

    sbuild = ensure_dir(work_dir / "build_solvent")
    s_prm = sbuild / "solvent.prmtop"; s_crd = sbuild / "solvent.inpcrd"
    s_pdb = sbuild / "solvent.pdb"; s_leap = sbuild / "tleap_solvent.in"
    s_leap.write_text(f"""# RBFE dual-topology SOLVENT leg (A->B).
source leaprc.{sys_cfg.ligand_ff}
source {water_leaprc}
loadAmberParams {frc_a.resolve()}
loadAmberParams {frc_b.resolve()}
{RESNAME_A} = loadMol2 {mol2_a.resolve()}
{RESNAME_B} = loadMol2 {mol2_b.resolve()}
system = combine {{ {RESNAME_A} {RESNAME_B} }}
{_solvate_cmd(sys_cfg, buf)}
{_ion_block(sys_cfg)}
saveAmberParm system {s_prm} {s_crd}
savePDB system {s_pdb}
quit
""")
    log.info("RBFE edge: building solvated dual-topology solvent leg (tleap)")
    _tleap_build(s_leap, sbuild, s_prm)

    out = {
        "complex_prmtop": str(c_prm), "complex_inpcrd": str(c_crd),
        "complex_pdb": str(c_pdb),
        "solvent_prmtop": str(s_prm), "solvent_inpcrd": str(s_crd),
        "timask1": f":{RESNAME_A}", "timask2": f":{RESNAME_B}",
        "scmask1": f":{RESNAME_A}", "scmask2": f":{RESNAME_B}",
        "ligand_a_file": str(ligand_a_file), "ligand_b_file": str(ligand_b_file),
    }
    (work_dir / "rbfe_edge_topology.json").write_text(json.dumps(out, indent=2))
    log.info("RBFE edge topology done: complex=%s solvent=%s", c_prm, s_prm)
    return out
