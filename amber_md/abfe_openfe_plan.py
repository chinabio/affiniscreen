"""amber_md/abfe_openfe_plan.py -- plan an ABSOLUTE binding free energy (ABFE)
calculation with OpenFE's AbsoluteBindingProtocol (runs on OpenMM).

Why this exists
---------------
The AMBER/pmemd ABFE path (amber_md/fep.py) requires hand-tuning GTI softcore
flags, scalpha/scbeta, barostat coupling, and lambda schedules -- and still
hits softcore singularities / box collapses. OpenFE's AbsoluteBindingProtocol
does ALL of that internally and is heavily validated:
  * Boresch orientational restraints (auto-selected) with lambda_restraints
  * charge annihilation + Lennard-Jones decoupling of the ligand
  * BOTH the complex and the solvent legs (the solvent leg is created by the
    protocol from the complex state -- you do NOT define it separately)
  * MBAR analysis via the protocol's gather()/get_estimate()

Mirrors the OpenFE `abfe_tutorial.ipynb` (toluene / T4 lysozyme) step-for-step.

Usage
-----
    python -m amber_md.abfe_openfe_plan \
        --protein protein.pdb --ligands ligand.sdf --out abfe_setup \
        --complex-ns 5.0 --solvent-ns 0.5 --platform CUDA \
        --charge-method am1bcc --repeats 1

One transformation JSON is written per ligand under
<out>/transformations/<ligand>.json, ready for `openfe quickrun`.

WARNING: OpenFE 1.x API attribute names shift between releases. Fragile settings
assignments are guarded; if one is missing on your build the planner logs a
warning and proceeds with the protocol default rather than crashing.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from openff.units import unit

import openfe
from openfe import (
    ProteinComponent, SolventComponent, SmallMoleculeComponent,
    ChemicalSystem, Transformation,
)
from openfe.protocols.openmm_afe import AbsoluteBindingProtocol


# --------------------------------------------------------------------------
# ligand loading + partial charge assignment (tutorial section 1-2)
# --------------------------------------------------------------------------
def _rdkit_mols_from_file(path: Path):
    """Yield RDKit mols from .sdf/.mol (SDMolSupplier) or .mol2 (per-block)."""
    from rdkit import Chem
    ext = path.suffix.lower()
    if ext == ".mol2":
        text = path.read_text(errors="ignore")
        marker = "@<TRIPOS>MOLECULE"
        if marker not in text:
            yield Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
            return
        for blk in (marker + b for b in text.split(marker)[1:]):
            m = None
            try:
                m = Chem.MolFromMol2Block(blk, removeHs=False, sanitize=True)
            except Exception:
                m = None
            if m is None:
                try:
                    from openbabel import pybel
                    ob = pybel.readstring("mol2", blk)
                    m = Chem.MolFromMolBlock(ob.write("mol"),
                                             removeHs=False, sanitize=True)
                except Exception:
                    m = None
            yield m
    else:
        for m in Chem.SDMolSupplier(str(path), removeHs=False):
            yield m


def load_ligands(sdf: Path):
    """Load ligands from .sdf / .mol / .mol2 into SmallMoleculeComponents."""
    out = []
    for i, m in enumerate(_rdkit_mols_from_file(sdf)):
        if m is None:
            print(f"  [warn] ligand record {i} failed to parse; skipping.",
                  file=sys.stderr)
            continue
        if not m.HasProp("_Name") or not m.GetProp("_Name").strip():
            m.SetProp("_Name", f"ligand_{i}")
        out.append(SmallMoleculeComponent.from_rdkit(m))
    if not out:
        raise SystemExit(f"No valid ligands parsed from {sdf}.")
    return out


def assign_charges(ligands, method="am1bcc", backend="ambertools"):
    """Pre-assign partial charges (tutorial section 2). Done once up-front so
    quickrun does not recompute am1bcc per leg. Best-effort: if the OpenFE
    charge utilities are unavailable on this build, return ligands unchanged
    and let the protocol charge them at runtime."""
    try:
        from openfe.protocols.openmm_utils.omm_settings import (
            OpenFFPartialChargeSettings)
        from openfe.protocols.openmm_utils.charge_generation import (
            bulk_assign_partial_charges)
    except Exception as e:  # pragma: no cover - version dependent
        print(f"  [warn] OpenFE charge utilities unavailable ({e}); "
              "protocol will assign charges at runtime.", file=sys.stderr)
        return ligands
    cs = OpenFFPartialChargeSettings(partial_charge_method=method,
                                     off_toolkit_backend=backend)
    charged = bulk_assign_partial_charges(
        ligands,
        overwrite=False,
        method=cs.partial_charge_method,
        toolkit_backend=cs.off_toolkit_backend,
        generate_n_conformers=getattr(cs, "number_of_conformers", None),
        nagl_model=getattr(cs, "nagl_model", None),
        processors=1,
    )
    return charged


# --------------------------------------------------------------------------
# protocol construction (tutorial section 4)
# --------------------------------------------------------------------------
def _try_set(obj_path, value, label):
    """Assign a possibly-renamed nested setting; warn (don't crash) if absent."""
    try:
        target, attr = obj_path
        setattr(target, attr, value)
    except Exception as e:  # pragma: no cover
        print(f"  [warn] could not set {label} ({e}); using protocol default.",
              file=sys.stderr)


def build_protocol(args) -> AbsoluteBindingProtocol:
    settings = AbsoluteBindingProtocol.default_settings()

    # one repeat per quickrun call so repeats parallelise as separate jobs
    settings.protocol_repeats = args.repeats

    # Boresch restraint host distances (tutorial: avoid PBC issues)
    _try_set((settings.restraint_settings, "host_min_distance"),
             args.host_min_nm * unit.nanometer, "host_min_distance")
    _try_set((settings.restraint_settings, "host_max_distance"),
             args.host_max_nm * unit.nanometer, "host_max_distance")

    # compute platform
    _try_set((settings.engine_settings, "compute_platform"),
             args.platform, "compute_platform")

    # small-molecule force field
    _try_set((settings.forcefield_settings, "small_molecule_forcefield"),
             args.forcefield, "small_molecule_forcefield")

    # production lengths for each leg (complex is the expensive one)
    _try_set((settings.complex_simulation_settings, "production_length"),
             args.complex_ns * unit.nanosecond, "complex production_length")
    _try_set((settings.solvent_simulation_settings, "production_length"),
             args.solvent_ns * unit.nanosecond, "solvent production_length")

    return AbsoluteBindingProtocol(settings=settings)


# --------------------------------------------------------------------------
# transformation building (tutorial section 3 + 4)
# --------------------------------------------------------------------------
def build_transformation(ligand, protein, solvent, protocol) -> Transformation:
    # state A: ligand fully interacting in the complex
    systemA = ChemicalSystem(
        {"ligand": ligand, "protein": protein, "solvent": solvent},
        name=ligand.name)
    # state B: ligand decoupled -> only protein + solvent.
    # The protocol derives the SOLVENT leg internally from these states.
    systemB = ChemicalSystem(
        {"protein": protein, "solvent": solvent})
    # ABFE has no atom mapping (nothing is mapped) -> mapping=None
    return Transformation(stateA=systemA, stateB=systemB, mapping=None,
                          protocol=protocol, name=ligand.name)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Plan OpenFE AbsoluteBindingProtocol (ABFE) transformations.")
    ap.add_argument("--protein", required=True, type=Path)
    ap.add_argument("--ligands", required=True, type=Path,
                    help="SDF/MOL/MOL2; one transformation is planned per record.")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--complex-ns", type=float, default=5.0,
                    help="Complex-leg production length (ns). Tutorial smoke "
                         "test uses 0.5; production default 10.")
    ap.add_argument("--solvent-ns", type=float, default=0.5,
                    help="Solvent-leg production length (ns).")
    ap.add_argument("--platform", default="CUDA", choices=["CUDA", "OpenCL", "CPU"])
    ap.add_argument("--forcefield", default="openff-2.2.0")
    ap.add_argument("--charge-method", default="am1bcc",
                    choices=["am1bcc", "am1bccelf10", "nagl", "espaloma"])
    ap.add_argument("--charge-backend", default="ambertools")
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--host-min-nm", type=float, default=0.5)
    ap.add_argument("--host-max-nm", type=float, default=1.5)
    args = ap.parse_args(argv)

    ligands = load_ligands(args.ligands)
    ligands = assign_charges(ligands, args.charge_method, args.charge_backend)
    protein = ProteinComponent.from_pdb_file(str(args.protein))
    solvent = SolventComponent()  # water + 0.15 M NaCl (OpenFE default)
    protocol = build_protocol(args)

    tdir = args.out / "transformations"
    tdir.mkdir(parents=True, exist_ok=True)

    n = 0
    for lig in ligands:
        tr = build_transformation(lig, protein, solvent, protocol)
        dest = tdir / f"{tr.name}.json"
        try:
            tr.dump(dest)
        except Exception:
            dest.write_text(tr.to_json())
        n += 1
        print(f"  planned ABFE transformation: {tr.name} -> {dest}")

    print(f"Planned {n} ABFE transformation(s) -> {tdir}")
    print("Run each with:  openfe quickrun <json> -o results.json -d workdir")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
