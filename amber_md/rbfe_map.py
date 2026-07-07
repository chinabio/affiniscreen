#!/usr/bin/env python3
"""rbfe_map.py - RBFE perturbation-map generator for Amber (two backends).

BACKENDS (auto; override --backend):
  * openfe : OpenFE validated planners (generate_minimal_redundant_network /
             generate_minimal_spanning_network / generate_radial_network)
             + LomapAtomMapper + lomap_scorers. Needs conda env with
             openfe+lomap+konnektor+gufe+rdkit.
  * rdkit  : dependency-light fallback (rdkit+networkx). PyAutoFEP-style MCS
             scoring (atoms hashed by ring/hybrid/Z; cost=1-exp(beta*pert),
             beta=ln(0.2)/median) + cycle-closure-aware redundant network.

Does NOT build prmtops; emits residue masks + per-edge atom mapping
(mapping.json) for a validated topology engine (tleap/pmx/femto), then a
fep_driver.py command scaffold (run_edges.sh).
License: derives from PyAutoFEP (GPLv2); references OpenFE (MIT) API -> GPLv2.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import argparse, csv, itertools, json, math, statistics, sys
from pathlib import Path
try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
except ImportError:
    sys.stderr.write("rbfe_map requires RDKit.\n"); raise
try:
    import networkx as nx
except ImportError:
    nx = None


def load_molecules(sdf_paths):
    molecules = {}
    for p in sdf_paths:
        for mol in Chem.SDMolSupplier(str(p), removeHs=False):
            if mol is None:
                continue
            name = (mol.GetProp("_Name")
                    if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                    else f"lig{len(molecules)+1}")
            mol.SetProp("_Name", name); molecules[name] = mol
    if len(molecules) < 2:
        raise RuntimeError("Need >=2 named molecules in input SDF(s).")
    return molecules


def suggest_masks(resname_a=":L1", resname_b=":L2"):
    return dict(timask1=resname_a, timask2=resname_b,
                scmask1=resname_a, scmask2=resname_b)


def openfe_available():
    try:
        import openfe  # noqa
        from openfe.setup.ligand_network_planning import (  # noqa
            generate_minimal_redundant_network)
        return True
    except Exception:
        return False


# ---------------- RDKit fallback ----------------
def _hash_atom(a):
    return 1000*int(a.IsInRing()) + 100*int(a.GetHybridization()) + a.GetAtomicNum()


def _clean_mcs_smarts(mol_a, hashed_smarts):
    patt = Chem.MolFromSmarts(hashed_smarts)
    if patt is None:
        return hashed_smarts
    a = Chem.RemoveHs(Chem.Mol(mol_a))
    match = a.GetSubstructMatch(patt)
    if not match:
        return hashed_smarts
    return Chem.MolFragmentToSmiles(a, atomsToUse=list(match), canonical=True)


def find_mcs(mol_a, mol_b, timeout=300):
    a = Chem.RemoveHs(Chem.Mol(mol_a)); b = Chem.RemoveHs(Chem.Mol(mol_b))
    if Chem.MolToSmiles(a) == Chem.MolToSmiles(b):
        smarts = Chem.MolToSmarts(a)
    else:
        ha, hb = Chem.Mol(a), Chem.Mol(b)
        for m in (ha, hb):
            for at in m.GetAtoms():
                at.SetIsotope(_hash_atom(at))
        res = rdFMCS.FindMCS([ha, hb],
                             atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             completeRingsOnly=True, ringMatchesRingOnly=True,
                             matchChiralTag=False, timeout=timeout)
        if res.canceled or res.numAtoms == 0:
            raise RuntimeError(f"MCS failed: {mol_a.GetProp('_Name')} vs "
                               f"{mol_b.GetProp('_Name')}")
        smarts = res.smartsString
    patt = Chem.MolFromSmarts(smarts)
    ma = a.GetSubstructMatch(patt); mb = b.GetSubstructMatch(patt)
    return smarts, patt.GetNumAtoms(), (list(zip(ma, mb)) if ma and mb else [])


def rdkit_cost_edges(molecules, use_hs=False, timeout=300, verbose=False):
    names = list(molecules); edges = []
    for ni, nj in itertools.combinations(names, 2):
        mi, mj = molecules[ni], molecules[nj]
        if mi.GetNumHeavyAtoms() < mj.GetNumHeavyAtoms():
            ni, nj, mi, mj = nj, ni, mj, mi
        smarts, core_heavy, mapping = find_mcs(mi, mj, timeout=timeout)
        if use_hs:
            ai, aj = mi.GetNumAtoms(), mj.GetNumAtoms()
            core = Chem.MolFromSmarts(smarts).GetNumAtoms()
        else:
            ai, aj = mi.GetNumHeavyAtoms(), mj.GetNumHeavyAtoms(); core = core_heavy
        perturbed = (ai-core)+(aj-core)
        if perturbed == 0:
            if verbose:
                sys.stderr.write(f"WARN: {ni}-{nj} perturbs 0 atoms; skip.\n")
            continue
        edges.append(dict(lig_a=ni, lig_b=nj, core_heavy=core_heavy,
                          perturbed_atoms=perturbed,
                          mcs_smarts=_clean_mcs_smarts(mi, smarts),
                          mapping=mapping, cost=None, score=None))
    if not edges:
        raise RuntimeError("No valid perturbation edges.")
    med = statistics.median(e["perturbed_atoms"] for e in edges)
    beta = math.log(0.2)/med if med else -1.0
    for e in edges:
        e["cost"] = 1.0 - math.exp(beta*e["perturbed_atoms"])
        e["score"] = 1.0 - e["cost"]
    return edges, beta, med


def _mst_only(edges, names):
    g = nx.Graph(); g.add_nodes_from(names)
    for e in edges:
        g.add_edge(e["lig_a"], e["lig_b"], weight=e["cost"], _e=e)
    if not nx.is_connected(g):
        raise RuntimeError("Ligand similarity graph is disconnected.")
    return [d["_e"] for _, _, d in nx.minimum_spanning_tree(g, weight="weight").edges(data=True)]


def _cycle_closure_network(edges, names, target_degree=2):
    g = nx.Graph(); g.add_nodes_from(names)
    for e in edges:
        g.add_edge(e["lig_a"], e["lig_b"], weight=e["cost"], _e=e)
    if not nx.is_connected(g):
        raise RuntimeError("Ligand similarity graph is disconnected.")
    mst = nx.minimum_spanning_tree(g, weight="weight")
    chosen = {frozenset((u, v)): d["_e"] for u, v, d in mst.edges(data=True)}
    deg = {n: 0 for n in names}
    for fs in chosen:
        for n in fs:
            deg[n] += 1
    for e in sorted((e for e in edges
                     if frozenset((e["lig_a"], e["lig_b"])) not in chosen),
                    key=lambda e: e["cost"]):
        a, b = e["lig_a"], e["lig_b"]
        if deg[a] < target_degree or deg[b] < target_degree:
            chosen[frozenset((a, b))] = e; deg[a] += 1; deg[b] += 1
    return list(chosen.values())


def _star_network(edges, names, hub=None):
    look = {frozenset((e["lig_a"], e["lig_b"])): e for e in edges}
    def tot(c):
        s = 0.0
        for o in names:
            if o == c:
                continue
            e = look.get(frozenset((c, o)))
            if e is None:
                return math.inf
            s += e["cost"]
        return s
    if hub is None:
        hub = min(names, key=tot)
        if tot(hub) == math.inf:
            raise RuntimeError("No node connects to all others; star impossible.")
    return [look[frozenset((hub, o))] for o in names if o != hub], hub


def run_rdkit_backend(molecules, map_type, hub, min_degree, use_hs, timeout, verbose):
    if nx is None:
        raise RuntimeError("networkx required for rdkit backend.")
    names = list(molecules)
    edges, beta, med = rdkit_cost_edges(molecules, use_hs, timeout, verbose)
    hub_out = None
    if map_type == "star":
        chosen, hub_out = _star_network(edges, names, hub)
    elif map_type == "mst":
        chosen = _mst_only(edges, names)
    else:
        chosen = _cycle_closure_network(edges, names, target_degree=min_degree)
    return chosen, dict(backend="rdkit", candidate_edges=len(edges),
                        median_perturbed=med, beta=beta, hub=hub_out)


# ---------------- OpenFE backend (lazy; validate in conda env) ----------------
def run_openfe_backend(molecules, map_type, hub, min_degree, verbose):
    from gufe import SmallMoleculeComponent
    from lomap import LomapAtomMapper
    from openfe.setup.atom_mapping import lomap_scorers
    from openfe.setup.ligand_network_planning import (
        generate_minimal_redundant_network, generate_minimal_spanning_network,
        generate_radial_network)
    try:
        from openfe.setup.ligand_network_planning import generate_lomap_network
    except Exception:
        generate_lomap_network = None

    comps = [SmallMoleculeComponent.from_rdkit(m) for m in molecules.values()]
    name_of = {id(c): n for c, n in zip(comps, molecules.keys())}
    mapper = LomapAtomMapper(); scorer = lomap_scorers.default_lomap_score

    if map_type == "star":
        if hub is None:
            raise RuntimeError("openfe star map needs --hub LIGNAME.")
        central = next(c for c, n in zip(comps, molecules) if n == hub)
        network = generate_radial_network(
            ligands=[c for c in comps if c is not central],
            central_ligand=central, mappers=[mapper], scorer=scorer)
    elif map_type == "mst":
        network = generate_minimal_spanning_network(ligands=comps, mappers=[mapper], scorer=scorer)
    elif map_type == "lomap" and generate_lomap_network is not None:
        network = generate_lomap_network(molecules=comps, mappers=[mapper], scorer=scorer)
    else:
        network = generate_minimal_redundant_network(
            ligands=comps, mappers=[mapper], scorer=scorer, mst_num=max(2, min_degree))

    chosen = []
    for edge in network.edges:
        a, b = edge.componentA, edge.componentB
        na = name_of.get(id(a), a.name); nb = name_of.get(id(b), b.name)
        try:
            sc = scorer(edge)
        except Exception:
            sc = None
        mapping = list(edge.componentA_to_componentB.items()) \
            if hasattr(edge, "componentA_to_componentB") else []
        chosen.append(dict(lig_a=na, lig_b=nb, core_heavy=len(mapping),
            perturbed_atoms=(a.to_rdkit().GetNumHeavyAtoms()
                             + b.to_rdkit().GetNumHeavyAtoms() - 2*len(mapping)),
            cost=(1.0-sc) if isinstance(sc, (int, float)) else None,
            score=sc, mcs_smarts="", mapping=mapping))
    return chosen, dict(backend="openfe", candidate_edges=None,
                        median_perturbed=None, beta=None, hub=hub)


# ---------------- output ----------------
def write_outputs(outdir, molecules, chosen, meta, resname_a, resname_b):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    masks = suggest_masks(resname_a, resname_b)
    with open(outdir/"nodes.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ligand", "heavy_atoms", "smiles"])
        for n, m in molecules.items():
            w.writerow([n, m.GetNumHeavyAtoms(), Chem.MolToSmiles(m)])
    with open(outdir/"edges.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge", "lig_a", "lig_b", "core_heavy", "perturbed_atoms",
                    "cost", "score", "timask1", "timask2", "scmask1", "scmask2", "mcs_smarts"])
        for e in chosen:
            w.writerow([f"{e['lig_a']}~{e['lig_b']}", e["lig_a"], e["lig_b"],
                        e.get("core_heavy", ""), e.get("perturbed_atoms", ""),
                        "" if e.get("cost") is None else f"{e['cost']:.4f}",
                        "" if e.get("score") is None else f"{e['score']:.4f}",
                        masks["timask1"], masks["timask2"], masks["scmask1"],
                        masks["scmask2"], e.get("mcs_smarts", "")])
    (outdir/"mapping.json").write_text(json.dumps(
        {f"{e['lig_a']}~{e['lig_b']}": [list(p) for p in e.get("mapping", [])]
         for e in chosen}, indent=2))
    if nx is not None:
        g = nx.Graph()
        for n, m in molecules.items():
            g.add_node(n, heavy_atoms=m.GetNumHeavyAtoms())
        for e in chosen:
            g.add_edge(e["lig_a"], e["lig_b"], perturbed_atoms=e.get("perturbed_atoms", 0),
                       cost=round(e["cost"], 4) if e.get("cost") is not None else 0.0)
        nx.write_graphml(g, str(outdir/"map.graphml"))
        cycles = g.number_of_edges()-g.number_of_nodes()+nx.number_connected_components(g)
        deg = dict(g.degree())
        (outdir/"diagnostics.txt").write_text(
            f"backend           : {meta.get('backend')}\n"
            f"nodes             : {g.number_of_nodes()}\n"
            f"edges             : {g.number_of_edges()}\n"
            f"independent cycles: {cycles}\n"
            f"min node degree   : {min(deg.values())}\n"
            f"max node degree   : {max(deg.values())}\n"
            f"isolated nodes    : {[n for n,d in deg.items() if d==0]}\n")
    with open(outdir/"run_edges.sh", "w") as f:
        f.write("#!/bin/bash\n# Auto-generated by rbfe_map.py\n")
        f.write("# Build a dual-topology prmtop per edge FIRST (tleap/pmx/femto), then:\n")
        for e in chosen:
            ed = f"{e['lig_a']}~{e['lig_b']}"
            f.write(f"python -m amber_md.fep_driver --mode legacy --work-dir ./edges/{ed} "
                    f"--complex-prmtop ./edges/{ed}/complex.prmtop "
                    f"--complex-inpcrd ./edges/{ed}/complex.inpcrd "
                    f"--solvent-prmtop ./edges/{ed}/solvent.prmtop "
                    f"--solvent-inpcrd ./edges/{ed}/solvent.inpcrd "
                    f"--timask1 '{masks['timask1']}' --timask2 '{masks['timask2']}' "
                    f"--scmask1 '{masks['scmask1']}' --scmask2 '{masks['scmask2']}' --submit\n")
    (outdir/"run_edges.sh").chmod(0o755)


def main(argv=None):
    p = argparse.ArgumentParser(prog="rbfe_map",
        description="RBFE perturbation-map generator for Amber (OpenFE + RDKit).")
    p.add_argument("-i", "--input", nargs="+", required=True)
    p.add_argument("-o", "--outdir", default="rbfe_map_out")
    p.add_argument("--backend", choices=["auto", "openfe", "rdkit"], default="auto")
    p.add_argument("--map-type", choices=["redundant", "mst", "star", "lomap"], default="redundant")
    p.add_argument("--hub", default=None)
    p.add_argument("--min-degree", type=int, default=2)
    p.add_argument("--use-hs", action="store_true")
    p.add_argument("--resname-a", default=":L1")
    p.add_argument("--resname-b", default=":L2")
    p.add_argument("--mcs-timeout", type=int, default=300)
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args(argv)
    molecules = load_molecules(a.input)
    backend = a.backend
    if backend == "auto":
        backend = "openfe" if openfe_available() else "rdkit"
    if backend == "openfe" and not openfe_available():
        sys.stderr.write("ERROR: openfe not importable; use --backend rdkit.\n"); return 2
    if backend == "rdkit" and a.map_type == "lomap":
        sys.stderr.write("ERROR: --map-type lomap requires --backend openfe.\n"); return 2
    if backend == "openfe":
        chosen, meta = run_openfe_backend(molecules, a.map_type, a.hub, a.min_degree, a.verbose)
    else:
        chosen, meta = run_rdkit_backend(molecules, a.map_type, a.hub, a.min_degree,
                                         a.use_hs, a.mcs_timeout, a.verbose)
    write_outputs(a.outdir, molecules, chosen, meta, a.resname_a, a.resname_b)
    print(f"backend         : {meta['backend']}")
    print(f"ligands         : {len(molecules)}")
    if meta.get("candidate_edges") is not None:
        print(f"candidate edges : {meta['candidate_edges']} "
              f"(median perturbed={meta['median_perturbed']}, beta={meta['beta']:.4f})")
    print(f"map type        : {a.map_type}" + (f"  hub={meta['hub']}" if meta.get('hub') else ""))
    print(f"chosen edges    : {len(chosen)}")
    print(f"outputs written : {Path(a.outdir).resolve()}")
    for e in sorted(chosen, key=lambda e: (e.get('cost') or 0)):
        c = "" if e.get("cost") is None else f"cost={e['cost']:.3f}"
        print(f"  {e['lig_a']:>10} ~ {e['lig_b']:<10}  perturbed={e.get('perturbed_atoms','?')!s:<3} {c}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
