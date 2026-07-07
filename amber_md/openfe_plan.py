"""amber_md/openfe_plan.py -- plan an RBFE network with custom SIMULATION
settings that the `openfe plan-rbfe-network` CLI cannot express.

The CLI `-s` YAML only covers mapper/network/charges. Simulation length,
compute platform, sampler, and small-molecule force field live in the protocol
object, so we build the network in Python.

Usage:
    python -m amber_md.openfe_plan \\
        --protein protein.pdb --ligands ligands.sdf --out network_setup \\
        --sim-ns 5.0 --equil-ns 1.0 --platform CUDA --sampler repex \\
        --forcefield openff-2.2.0 --network mst --charges am1bcc

Each transformation is written to <out>/transformations/<edge>.json, ready for
`openfe quickrun`. Planned with n_repeats=1 so repeats parallelise.

WARNING: The OpenFE Python API import paths and attribute names shift between
1.x releases. On a dev build, some attributes below may differ. Fragile bits are
guarded with try/except. Test once and adjust the marked lines if needed; for a
first smoke test the plain CLI planner is the safest path.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations

import argparse
from pathlib import Path

from openff.units import unit

from openfe import (
    ProteinComponent, SolventComponent, SmallMoleculeComponent,
    ChemicalSystem, Transformation, AlchemicalNetwork,
)
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.setup import KartografAtomMapper
from openfe.setup.ligand_network_planning import (
    generate_minimal_spanning_network,
    generate_minimal_redundant_network,
    generate_radial_network,
)

try:
    from kartograf.atom_mapping_scorer import default_lomap_score as _SCORE
except Exception:  # pragma: no cover - version dependent
    _SCORE = None

from rdkit import Chem


def load_ligands(sdf: Path):
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    out = []
    for i, m in enumerate(suppl):
        if m is None:
            continue
        if not m.HasProp("_Name") or not m.GetProp("_Name").strip():
            m.SetProp("_Name", f"ligand_{i}")
        out.append(SmallMoleculeComponent.from_rdkit(m))
    if len(out) < 2:
        raise SystemExit(f"Need >=2 valid ligands for RBFE; got {len(out)}.")
    return out


def build_protocol(args) -> RelativeHybridTopologyProtocol:
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.protocol_repeats = 1  # one repeat per quickrun call
    settings.simulation_settings.equilibration_length = args.equil_ns * unit.nanosecond
    settings.simulation_settings.production_length = args.sim_ns * unit.nanosecond
    try:
        settings.engine_settings.compute_platform = args.platform
    except Exception:
        pass
    try:
        settings.forcefield_settings.small_molecule_forcefield = args.forcefield
    except Exception:
        pass
    try:
        settings.simulation_settings.sampler_method = args.sampler
    except Exception:
        pass
    return RelativeHybridTopologyProtocol(settings)


def plan_network(ligs, mapper, network_algo):
    kw = {"mappers": [mapper]}
    if _SCORE is not None:
        kw["scorer"] = _SCORE
    if network_algo == "mst":
        return generate_minimal_spanning_network(ligs, **kw)
    if network_algo == "redundant":
        return generate_minimal_redundant_network(ligs, **kw)
    if network_algo == "radial":
        return generate_radial_network(ligs, central_ligand=ligs[0], **kw)
    raise SystemExit(f"Unknown network algo: {network_algo}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True, type=Path)
    ap.add_argument("--ligands", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sim-ns", type=float, default=5.0)
    ap.add_argument("--equil-ns", type=float, default=1.0)
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--sampler", default="repex")
    ap.add_argument("--forcefield", default="openff-2.2.0")
    ap.add_argument("--network", default="mst",
                    choices=["mst", "redundant", "radial"])
    ap.add_argument("--charges", default="am1bcc")  # informational here
    args = ap.parse_args()

    ligs = load_ligands(args.ligands)
    protein = ProteinComponent.from_pdb_file(str(args.protein))
    solvent = SolventComponent()  # water + 0.15 M NaCl (OpenFE default)
    protocol = build_protocol(args)
    mapper = KartografAtomMapper()

    ligand_network = plan_network(ligs, mapper, args.network)

    transformations = []
    for mapping in ligand_network.edges:
        a, b = mapping.componentA, mapping.componentB
        sysA_c = ChemicalSystem({"ligand": a, "protein": protein, "solvent": solvent})
        sysB_c = ChemicalSystem({"ligand": b, "protein": protein, "solvent": solvent})
        sysA_s = ChemicalSystem({"ligand": a, "solvent": solvent})
        sysB_s = ChemicalSystem({"ligand": b, "solvent": solvent})
        for leg, (sA, sB) in {"complex": (sysA_c, sysB_c),
                              "solvent": (sysA_s, sysB_s)}.items():
            transformations.append(Transformation(
                stateA=sA, stateB=sB, protocol=protocol,
                mapping={"ligand": mapping},
                name=f"{a.name}_{b.name}_{leg}"))

    network = AlchemicalNetwork(transformations)
    out = args.out
    (out / "transformations").mkdir(parents=True, exist_ok=True)
    try:
        network.to_json(out / "alchemical_network.json")
    except Exception:
        pass
    for t in network.edges:
        try:
            t.dump(out / "transformations" / f"{t.name}.json")
        except Exception:
            (out / "transformations" / f"{t.name}.json").write_text(t.to_json())
    print(f"Planned {len(transformations)} transformations "
          f"({len(ligand_network.edges)} edges x 2 legs) -> {out}")


if __name__ == "__main__":
    main()
