"""Shared widgets and helpers for the ABFE and RBFE Streamlit pages.

Lives at amber_md/gui/fep_common.py.

Both pages import from here so the protein/ligand input widget, MD-settings
widget, HPC-settings widget, and SDF parsing live in exactly one place.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import streamlit as st

from amber_md.gui.common import file_picker


# ---------------------------------------------------------------------------
# Data classes returned by the widgets
# ---------------------------------------------------------------------------
@dataclass
class CommonInputs:
    work_dir: Path
    protein_pdb: Path | None
    ligand_sdf: Path | None
    ligand_resname: str


@dataclass
class MDSettings:
    lambdas: list[float]
    nstlim_eq: int
    nstlim_prod: int
    dt: float
    cutoff: float
    temperature: float
    ntpr: int
    ntwx: int

    def as_cli(self) -> list[str]:
        return [
            "--lambdas", *[f"{l:.3f}" for l in self.lambdas],
            "--nstlim-eq",   str(self.nstlim_eq),
            "--nstlim-prod", str(self.nstlim_prod),
            "--dt",          str(self.dt),
            "--cutoff",      str(self.cutoff),
            "--ntpr",        str(self.ntpr),
            "--ntwx",        str(self.ntwx),
            "--temperature", str(self.temperature),
        ]


@dataclass
class HPCSettings:
    project: str
    queue: str
    walltime: str
    n_gpu: int
    modules: list[str]
    venv: str = ""
    use_hremd: bool = False
    exchange_freq: int = 1000

    def as_cli(self) -> list[str]:
        out = [
            "--project", self.project,
            "--queue", self.queue,
            "--walltime", self.walltime,
            "--n-gpu", str(self.n_gpu),
            "--modules", *self.modules,
        ]
        if self.venv.strip():
            out += ["--venv", self.venv.strip()]
        if self.use_hremd:
            out += ["--hremd", "--exchange-freq", str(self.exchange_freq)]
        return out


# ---------------------------------------------------------------------------
# SDF parsing
# ---------------------------------------------------------------------------
@dataclass
class LigandRecord:
    index: int          # 0-based position in the SDF
    name: str           # _Name property or fallback "ligand_<idx>"
    smiles: str | None  # may be None if RDKit fails
    n_heavy: int        # heavy-atom count, used as a quick sanity check


def _records_from_mols(mols) -> list[LigandRecord]:
    """Turn an iterable of (possibly-None) RDKit mols into LigandRecords."""
    from rdkit import Chem
    out: list[LigandRecord] = []
    for i, mol in enumerate(mols):
        if mol is None:
            out.append(LigandRecord(i, f"ligand_{i}_INVALID", None, 0))
            continue
        name = mol.GetProp("_Name") if (mol.HasProp("_Name")
                                        and mol.GetProp("_Name").strip()) else f"ligand_{i}"
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            smi = None
        out.append(LigandRecord(i, name, smi, mol.GetNumHeavyAtoms()))
    return out


def _load_mol2_records(path: Path) -> list[LigandRecord]:
    """Load one-or-many ligands from a MOL2 file.

    RDKit's Mol2 reader handles ONE molecule per call and is fragile, so we
    split the file on '@<TRIPOS>MOLECULE' first and parse each block. Falls
    back to OpenBabel via pybel if RDKit rejects a block.
    """
    from rdkit import Chem
    text = path.read_text(errors="ignore")
    marker = "@<TRIPOS>MOLECULE"
    if marker not in text:
        # not a real mol2; let RDKit try the whole file once
        m = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
        return _records_from_mols([m])
    blocks = [marker + b for b in text.split(marker)[1:]]
    mols = []
    for blk in blocks:
        m = None
        try:
            m = Chem.MolFromMol2Block(blk, removeHs=False, sanitize=True)
        except Exception:
            m = None
        if m is None:
            # OpenBabel fallback (handles SYBYL atom types RDKit chokes on)
            try:
                from openbabel import pybel
                ob = pybel.readstring("mol2", blk)
                m = Chem.MolFromMolBlock(ob.write("mol"), removeHs=False,
                                         sanitize=True)
            except Exception:
                m = None
        mols.append(m)
    return _records_from_mols(mols)


def load_sdf_records(sdf_path: Path) -> list[LigandRecord]:
    """Parse a ligand file (.sdf / .mol / .mol2) into LigandRecord objects.

    Despite the historical name, this now dispatches on file extension:
      * .sdf / .mol -> Chem.SDMolSupplier (multi-record)
      * .mol2       -> per-block MOL2 parse (RDKit, OpenBabel fallback)
    Tries RDKit first; falls back to a primitive record count so the page
    still loads (with a warning) on machines without RDKit.
    """
    if sdf_path is None or not sdf_path.exists():
        return []

    ext = sdf_path.suffix.lower()
    try:
        from rdkit import Chem
        if ext == ".mol2":
            return _load_mol2_records(sdf_path)
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
        return _records_from_mols(suppl)
    except ImportError:
        text = sdf_path.read_text(errors="ignore")
        if ext == ".mol2":
            n = text.count("@<TRIPOS>MOLECULE") or 1
        else:
            n = text.count("$$$$\n") or text.count("$$$$\r\n") or 1
        st.warning("RDKit not available; ligand metadata is limited.")
        return [LigandRecord(i, f"ligand_{i}", None, 0) for i in range(n)]


# ---------------------------------------------------------------------------
# Widgets
# ---------------------------------------------------------------------------
def render_common_inputs(key_prefix: str,
                         default_wd: str = "~/Run_dir/fep_run",
                         default_resname: str = "LIG") -> CommonInputs:
    """Sidebar widget: protein PDB + ligand SDF + work dir + resname."""
    with st.sidebar:
        st.header("Inputs")
        wd_str = st.text_input("Work directory", default_wd,
                               key=f"{key_prefix}_wd")
        wd = Path(wd_str).expanduser()

        protein_pdb = file_picker(
            "Protein PDB", f"{key_prefix}_protein",
            {".pdb"}, default_dir="~/Run_dir")
        ligand_sdf = file_picker(
            "Ligand SDF (1 or many records)", f"{key_prefix}_ligands",
            {".sdf", ".mol", ".mol2"}, default_dir="~/Run_dir")
        lig_resname = st.text_input(
            "Ligand residue name (used in masks)",
            default_resname, key=f"{key_prefix}_resname")

    return CommonInputs(
        work_dir=wd,
        protein_pdb=Path(protein_pdb) if protein_pdb else None,
        ligand_sdf=Path(ligand_sdf) if ligand_sdf else None,
        ligand_resname=lig_resname,
    )


def render_md_settings(key_prefix: str) -> MDSettings:
    """Tab body: lambda schedule + integrator + thermostat."""
    n_lam = st.slider("Number of lambda windows", 5, 24, 11,
                      key=f"{key_prefix}_nlam")
    spacing = st.selectbox("Spacing", ["linear", "denser at endpoints"],
                           key=f"{key_prefix}_spacing")
    if spacing == "linear":
        lambdas = [round(i / (n_lam - 1), 3) for i in range(n_lam)]
    else:
        lambdas = [round(0.5 * (1 - math.cos(math.pi * i / (n_lam - 1))), 3)
                   for i in range(n_lam)]
    st.code(", ".join(f"{l:.3f}" for l in lambdas))

    c1, c2, c3 = st.columns(3)
    nstlim_eq   = c1.number_input("nstlim (eq)",   1000, 5_000_000,   250_000, 1000,
                                  key=f"{key_prefix}_eq")
    nstlim_prod = c2.number_input("nstlim (prod)", 1000, 50_000_000, 20_000_000, 1000,
                                  key=f"{key_prefix}_prod")
    dt          = c3.number_input("dt (ps)", 0.0005, 0.004, 0.001, 0.0005,
                                  format="%.4f", key=f"{key_prefix}_dt")
    c1, c2, c3 = st.columns(3)
    cutoff      = c1.number_input("Cutoff (A)", 6.0, 14.0, 10.0, 0.5,
                                  key=f"{key_prefix}_cutoff")
    temperature = c2.number_input("Temperature (K)", 250.0, 350.0, 298.0, 1.0,
                                  key=f"{key_prefix}_T")
    ntpr        = c3.number_input("ntpr", 100, 100_000, 1000, 100,
                                  key=f"{key_prefix}_ntpr")
    ntwx = st.number_input("ntwx", 100, 100_000, 5000, 100,
                           key=f"{key_prefix}_ntwx")

    return MDSettings(lambdas=lambdas,
                      nstlim_eq=int(nstlim_eq), nstlim_prod=int(nstlim_prod),
                      dt=float(dt), cutoff=float(cutoff),
                      temperature=float(temperature),
                      ntpr=int(ntpr), ntwx=int(ntwx))


def render_hpc_settings(key_prefix: str,
                        allow_hremd: bool = True,
                        n_lambdas: int = 11) -> HPCSettings:
    """Tab body: LSF + HREMD knobs."""
    c1, c2 = st.columns(2)
    project  = c1.text_input("LSF project", "your-project", key=f"{key_prefix}_proj")
    queue    = c2.text_input("GPU queue", "gpu", key=f"{key_prefix}_q")
    walltime = c1.text_input("Walltime", "24:00", key=f"{key_prefix}_wt")
    n_gpu    = c2.number_input("GPUs per window (non-HREMD)", 1, 8, 1,
                               key=f"{key_prefix}_ngpu")
    modules  = st.text_input("module load (space-separated)", "amber/22",
                             key=f"{key_prefix}_mods").split()
    venv     = st.text_input(
        "venv activate script (optional)", "",
        placeholder="~/envs/affiniscreen/activate_amber_md.sh",
        help="Optional path to an activate script sourced before each job. "
             "Leave blank to use the current environment. '~' is expanded.",
        key=f"{key_prefix}_venv")

    use_hremd, exch_freq = False, 1000
    if allow_hremd:
        st.divider()
        st.markdown("### Hamiltonian replica exchange")
        use_hremd = st.checkbox(
            "Enable HREMD (one GPU per lambda)", value=False,
            key=f"{key_prefix}_hremd")
        exch_freq = st.number_input(
            "Exchange attempt period (steps)", 100, 50_000, 1000, 100,
            disabled=not use_hremd, key=f"{key_prefix}_xfreq")
        if use_hremd:
            st.info(f"HREMD will request **{n_lambdas}** GPUs per job.")

    return HPCSettings(project=project, queue=queue, walltime=walltime,
                       n_gpu=int(n_gpu), modules=modules, venv=venv,
                       use_hremd=use_hremd, exchange_freq=int(exch_freq))


# ---------------------------------------------------------------------------
# FEP-map planning (RBFE only)
# ---------------------------------------------------------------------------
@dataclass
class FEPEdge:
    a_idx: int
    b_idx: int
    a_name: str
    b_name: str
    score: float | None = None   # MCS-based, higher = more similar

@dataclass
class FEPMap:
    ligands: list[LigandRecord]
    edges: list[FEPEdge]

    @property
    def n_nodes(self) -> int:
        return len(self.ligands)

    @property
    def n_edges(self) -> int:
        return len(self.edges)


def plan_fep_map(ligands: list[LigandRecord],
                 sdf_path: Path,
                 strategy: str = "lomap") -> FEPMap:
    """Plan a perturbation network for RBFE.

    Tries LOMAP if installed (industry standard). Falls back to a
    nearest-neighbour-by-heavy-atom-count chain so the page is usable
    without LOMAP — clearly not as good, but never silently broken.
    """
    if len(ligands) < 2:
        return FEPMap(ligands=ligands, edges=[])

    try:
        import lomap  # type: ignore
        db = lomap.DBMolecules(str(sdf_path.parent),
                               output=False, parallel=1, verbose="off",
                               time=20, max3d=1.0, threed=False)
        # Restrict to our SDF
        # (LOMAP scans a directory; for a single SDF the simple path is to
        # pass the directory containing it and rely on it being the only SDF.)
        strict, loose, _ = db.build_matrices()
        nx_graph = db.build_graph()
        edges = []
        for u, v, data in nx_graph.edges(data=True):
            score = float(data.get("similarity", 0.0))
            edges.append(FEPEdge(
                a_idx=int(u), b_idx=int(v),
                a_name=ligands[int(u)].name, b_name=ligands[int(v)].name,
                score=score))
        return FEPMap(ligands=ligands, edges=edges)
    except ImportError:
        st.warning(
            "LOMAP is not installed; falling back to a heavy-atom-count "
            "nearest-neighbour chain. Install `lomap2` for production maps."
        )
        order = sorted(range(len(ligands)), key=lambda i: ligands[i].n_heavy)
        edges = [
            FEPEdge(a_idx=order[i], b_idx=order[i + 1],
                    a_name=ligands[order[i]].name,
                    b_name=ligands[order[i + 1]].name,
                    score=None)
            for i in range(len(order) - 1)
        ]
        # Add a couple of "long" edges so we get at least one cycle for QC.
        if len(order) >= 4:
            edges.append(FEPEdge(order[0], order[-1],
                                 ligands[order[0]].name,
                                 ligands[order[-1]].name, None))
            edges.append(FEPEdge(order[1], order[-2],
                                 ligands[order[1]].name,
                                 ligands[order[-2]].name, None))
        return FEPMap(ligands=ligands, edges=edges)


# ---------------------------------------------------------------------------
# Network solver (RBFE only)
# ---------------------------------------------------------------------------
def solve_network(edges_with_ddg: list[tuple[str, str, float, float]],
                  reference: str | None = None) -> dict[str, tuple[float, float]]:
    """Solve for per-ligand DeltaG from per-edge DDeltaG via MLE.

    Tries `cinnabar.FEMap` first; falls back to a plain weighted-least-squares
    implementation so the page works without cinnabar.

    Parameters
    ----------
    edges_with_ddg : list of (a_name, b_name, ddG_b_minus_a, err)
    reference : ligand name to pin to 0.0 (default: first ligand seen)
    """
    if not edges_with_ddg:
        return {}

    try:
        from cinnabar import FEMap, Measurement  # type: ignore
        fem = FEMap()
        for a, b, ddg, err in edges_with_ddg:
            fem.add_measurement(Measurement(
                labelA=a, labelB=b, DG=ddg, uncertainty=err,
                computational=True))
        fem.generate_absolute_values()
        out = {}
        for node in fem.absolute_dataframe.itertuples():
            out[node.label] = (float(node.DG), float(node.uncertainty))
        return out
    except ImportError:
        # ----- minimal MLE fallback -----
        import numpy as np
        names = sorted({n for a, b, _, _ in edges_with_ddg for n in (a, b)})
        idx = {n: i for i, n in enumerate(names)}
        ref = reference or names[0]
        n = len(names)
        m = len(edges_with_ddg)
        A = np.zeros((m + 1, n))
        y = np.zeros(m + 1)
        w = np.ones(m + 1)
        for k, (a, b, ddg, err) in enumerate(edges_with_ddg):
            A[k, idx[b]] = +1
            A[k, idx[a]] = -1
            y[k] = ddg
            w[k] = 1.0 / max(err, 1e-3) ** 2
        # pin the reference
        A[m, idx[ref]] = 1.0
        y[m] = 0.0
        w[m] = 1e6
        W = np.diag(w)
        AtA = A.T @ W @ A
        Aty = A.T @ W @ y
        x = np.linalg.solve(AtA, Aty)
        cov = np.linalg.inv(AtA)
        return {nm: (float(x[i]), float(np.sqrt(cov[i, i])))
                for nm, i in idx.items()}


def cycle_closure_residuals(edges_with_ddg, dg_per_ligand) -> list[tuple[list[str], float]]:
    """Find independent cycles and report the residual sum of DDeltaG.

    Returns list of (cycle_node_names, residual_kcal_mol). Residual should be
    near zero for a well-converged map; |residual| > 1 kcal/mol on any cycle
    is a red flag for that subgraph.
    """
    try:
        import networkx as nx
        G = nx.Graph()
        for a, b, ddg, err in edges_with_ddg:
            G.add_edge(a, b, ddg=ddg, err=err)
        out = []
        for cycle in nx.cycle_basis(G):
            residual = 0.0
            for u, v in zip(cycle, cycle[1:] + cycle[:1]):
                data = G[u][v]
                # signed: traverse u -> v in the direction it was added
                # we don't track direction here, so use |ddg|; this is a
                # rough QC, not a publication-grade closure
                residual += data["ddg"] if cycle.index(u) < cycle.index(v) else -data["ddg"]
            out.append((cycle, residual))
        return out
    except ImportError:
        return []
