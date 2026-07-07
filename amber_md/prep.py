"""PDB cleaning + ligand parametrization (v2.4.6).

v2.4.6 changes:
  * LigandParametrizer rejects multi-record SDF/MOL/MOL2 files BEFORE
    calling antechamber, with an actionable error message. Previously
    antechamber would fail deep inside with a misleading
    'Invalid number of atoms (0); MDL SDF supports at most 999 atoms.'
    which made batch-screen multi-SDFs accidentally fed to run_amber.py
    look like data corruption.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .utils import run, which_or_die, ensure_dir, CommandError
from .logger import get_logger
from .config import SystemConfig
log = get_logger()


# ---------------------------------------------------------------------
# Format detection (used by both PDBCleaner and LigandParametrizer)
# ---------------------------------------------------------------------
LIGAND_FORMATS = {
    ".sdf":  "sdf",
    ".mol":  "mdl",
    ".mol2": "mol2",
    ".pdb":  "pdb",
    ".xyz":  "xyz",
}

PROTEIN_FORMATS = {".pdb", ".mol2"}


def detect_ligand_format(path: Path) -> str:
    """Return the antechamber -fi value for a ligand file."""
    ext = Path(path).suffix.lower()
    if ext not in LIGAND_FORMATS:
        raise ValueError(
            f"Unsupported ligand format {ext!r}. "
            f"Supported: {sorted(LIGAND_FORMATS)}")
    return LIGAND_FORMATS[ext]


# ---------------------------------------------------------------------
# v2.4.6: multi-record detection for SDF/MOL/MOL2
# ---------------------------------------------------------------------
def _count_sdf_records(path: Path) -> int:
    """Count SDF/MOL records by counting '$$$$' delimiter lines.

    Per the MDL Molfile / SDF spec, each molecule ends with a line that
    is exactly '$$$$'. Empty SDFs (no records) return 0; valid single-
    molecule SDFs return 1; multi-molecule (concatenated) SDFs return >1.
    """
    n = 0
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                if line.strip() == "$$$$":
                    n += 1
    except OSError:
        return -1   # unreadable; let antechamber report it
    return n


def _count_mol2_records(path: Path) -> int:
    """Count MOL2 records by counting '@<TRIPOS>MOLECULE' headers."""
    n = 0
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                if line.startswith("@<TRIPOS>MOLECULE"):
                    n += 1
    except OSError:
        return -1
    return n


def _guard_single_record_ligand(ligand_file: Path) -> None:
    """Raise CommandError with actionable text if the ligand file holds
    more than one molecule. Silently accepts:
      - any format we can't sensibly count (PDB, XYZ)
      - files where the count is ambiguous (returns -1 on read errors)
      - files with exactly 1 record
      - SDF files with 0 '$$$$' markers (some writers omit the trailing
        delimiter on the last record; antechamber tolerates this case)
    """
    ext = ligand_file.suffix.lower()

    if ext in (".sdf", ".mol"):
        n = _count_sdf_records(ligand_file)
        if n > 1:
            raise CommandError(
                f"{ligand_file} contains {n} molecules.\n"
                f"run_amber.py / LigandParametrizer accepts a SINGLE-ligand file.\n"
                f"\n"
                f"Fixes:\n"
                f"  * For batch screening: use the Batch Screen page "
                f"(or `python -m amber_md.batch`) which auto-splits the SDF.\n"
                f"  * To run one ligand from this multi-SDF, extract the first record:\n"
                f"      awk '/^\\$\\$\\$\\$/{{print; exit}} {{print}}' "
                f"{ligand_file} > one_ligand.sdf\n"
                f"  * Or split all records with RDKit / obabel into per-molecule files.")
        return

    if ext == ".mol2":
        n = _count_mol2_records(ligand_file)
        if n > 1:
            raise CommandError(
                f"{ligand_file} contains {n} MOL2 molecules.\n"
                f"run_amber.py / LigandParametrizer accepts a SINGLE-ligand file.\n"
                f"\n"
                f"Fixes:\n"
                f"  * For batch screening: use the Batch Screen page "
                f"(or `python -m amber_md.batch`).\n"
                f"  * Split with obabel:\n"
                f"      obabel {ligand_file} -O lig_.mol2 --separate")
        return

    # .pdb / .xyz: a multi-model PDB or trajectory-style XYZ is a different
    # failure mode that antechamber will report clearly; not worth guarding.
    return


# ---------------------------------------------------------------------
# PDBCleaner — direct-input mode (skips the combined split)
# ---------------------------------------------------------------------
class PDBCleaner:
    def __init__(self, work_dir):
        self.work_dir = ensure_dir(work_dir)

    def clean(self, in_pdb, out_pdb=None):
        which_or_die("pdb4amber")
        in_pdb = Path(in_pdb).resolve()
        out_pdb = Path(out_pdb) if out_pdb else self.work_dir/f"{in_pdb.stem}_clean.pdb"
        run(["pdb4amber", "-i", str(in_pdb), "-o", str(out_pdb), "--reduce", "--dry"],
            cwd=self.work_dir)
        return out_pdb

    def clean_protein_only(self, in_protein, out_pdb=None):
        """Clean a protein-only file (no ligand expected).

        Accepts .pdb (passed directly to pdb4amber) or .mol2 (converted
        via antechamber first). The result is always a clean .pdb suitable
        for tleap's loadPDB.
        """
        in_protein = Path(in_protein).resolve()
        ext = in_protein.suffix.lower()
        if ext == ".pdb":
            return self.clean(in_protein, out_pdb)
        if ext == ".mol2":
            which_or_die("antechamber")
            intermediate = self.work_dir/f"{in_protein.stem}_from_mol2.pdb"
            run(["antechamber", "-i", str(in_protein), "-fi", "mol2",
                 "-o", str(intermediate), "-fo", "pdb", "-pf", "y"],
                cwd=self.work_dir)
            return self.clean(intermediate, out_pdb)
        raise ValueError(
            f"Unsupported protein format {ext!r}. Use .pdb or .mol2.")

    @staticmethod
    def split_complex(clean_pdb, ligand_resname, out_dir):
        """Split a combined complex PDB into protein.pdb + ligand.pdb.

        Bug fix (v2.2.2): END writes were OUTSIDE the with-block.
        Improved error message (v2.2.2): hints at the most common cause.
        """
        out_dir = ensure_dir(out_dir)
        prot = out_dir/"protein.pdb"
        lig  = out_dir/"ligand.pdb"
        resname_u = ligand_resname.upper()

        with open(clean_pdb) as f, open(prot, "w") as fp, open(lig, "w") as fl:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    if line[17:20].strip().upper() == resname_u:
                        fl.write(line)
                    else:
                        fp.write(line)
                elif line.startswith(("TER", "END")):
                    fp.write(line)
            fl.write("END\n")
            fp.write("END\n")

        if prot.stat().st_size <= 5:
            raise RuntimeError(
                f"split_complex: protein.pdb is empty. "
                f"Check that {clean_pdb} contains ATOM records.")
        if lig.stat().st_size <= 5:
            found = set()
            with open(clean_pdb) as f:
                for line in f:
                    if line.startswith("HETATM"):
                        found.add(line[17:20].strip())
            hint = (f" HETATM resnames present: {sorted(found)}. "
                    f"Use --lig-resname <one of those>." if found
                    else " No HETATM records present.")
            raise RuntimeError(
                f"split_complex: ligand.pdb is empty. "
                f"No residue named '{ligand_resname}' found in {clean_pdb}." + hint +
                f" Or supply protein and ligand as separate files: "
                f"--protein-file <prot.pdb> --ligand-file <lig.sdf|.mol2|.pdb>.")
        return prot, lig


# ---------------------------------------------------------------------
# LigandParametrizer — format-aware, multi-record-aware (v2.4.6)
# ---------------------------------------------------------------------
class LigandParametrizer:
    def __init__(self, work_dir, cfg):
        self.work_dir = ensure_dir(work_dir/"ligand")
        self.cfg = cfg

    def parametrize(self, ligand_file, resname="LIG"):
        """Run antechamber + parmchk2 on a SINGLE-ligand file.

        Accepts .pdb, .sdf, .mol, .mol2, .xyz — auto-detected from extension.
        For SDF / MOL / MOL2, rejects multi-record files up-front with an
        actionable error (v2.4.6).

        Charge method, net charge, multiplicity come from SystemConfig.
        """
        which_or_die("antechamber"); which_or_die("parmchk2")
        ligand_file = Path(ligand_file).resolve()
        in_fmt = detect_ligand_format(ligand_file)

        # v2.4.6: catch multi-record SDF/MOL/MOL2 before antechamber gets confused.
        _guard_single_record_ligand(ligand_file)

        mol2   = self.work_dir/f"{resname}.mol2"
        frcmod = self.work_dir/f"{resname}.frcmod"
        atype  = "gaff2" if self.cfg.ligand_ff.lower() == "gaff2" else "gaff"

        log.info("Parametrizing ligand %s (format=%s, resname=%s, "
                 "charge=%d, mult=%d, method=%s, atom_type=%s)",
                 ligand_file.name, in_fmt, resname,
                 self.cfg.ligand_charge, self.cfg.ligand_multiplicity,
                 self.cfg.charge_method, atype)

        run(["antechamber",
             "-i", str(ligand_file), "-fi", in_fmt,
             "-o", str(mol2),        "-fo", "mol2",
             "-c", self.cfg.charge_method,
             "-nc", str(self.cfg.ligand_charge),
             "-m",  str(self.cfg.ligand_multiplicity),
             "-at", atype, "-rn", resname,
             "-pf", "y", "-s", "2"], cwd=self.work_dir)
        run(["parmchk2", "-i", str(mol2), "-f", "mol2",
             "-o", str(frcmod), "-s", atype], cwd=self.work_dir)
        return mol2, frcmod