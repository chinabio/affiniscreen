"""
mmgbsa_openmm.py  --  OpenMM-MD + AmberTools MMPBSA.py MM-GBSA (Interpretation A).

Generates the production trajectory with OpenMM (free/open-source), then scores
binding free energy with AmberTools `MMPBSA.py` -- reusing the project's existing
TopologySplitter (amber_md.topology) and MMGBSAAnalyzer (amber_md.mmgbsa). No
dependence on the paid Amber `pmemd` engine.

Pipeline:
  1. Parametrize protein+ligand -> a single SOLVATED complex prmtop/inpcrd via
     tleap (AmberTools: pdb4amber, antechamber, parmchk2, tleap).
  2. OpenMM MD: createSystem from the Amber prmtop, minimize/equilibrate/produce,
     write a NetCDF trajectory whose atom order MATCHES solvated.prmtop.
  3. Score: TopologySplitter.split(solvated.prmtop) -> complex/receptor/ligand
     topologies; MMGBSAAnalyzer.run() -> FINAL_RESULTS_MMPBSA.dat (+ HTML report).

Design choices vs. the reference Colab recipe (quantaosun/Ambertools-OpenMM-MD):
  * Modern `import openmm` namespace (not deprecated `simtk.openmm`).
  * subprocess.run with explicit error checking (no bare os.system).
  * A sanitized AmberTools environment is the caller's responsibility (the GUI
    passes one); we do NOT source amber.sh or mutate AMBERHOME here.
  * Trajectory written as NetCDF (Amber-native) so MMPBSA.py -y reads it directly
    -- no cpptraj DCD->mdcrd conversion step.

CLI:
  python -m amber_md.mmgbsa_openmm \\
      --protein protein.pdb --ligand ligand.sdf --workdir run/ \\
      --prod-ns 10 --equil-ns 1 --igb 8 --salt 0.15 --platform CUDA \\
      [--ligand-charge 0] [--ligand-resname LIG] [--dry-run]
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from .logger import get_logger

log = get_logger()

# Tools we depend on (AmberTools). Checked in preflight().
REQUIRED_TOOLS = ["pdb4amber", "antechamber", "parmchk2", "tleap", "MMPBSA.py"]


# ---------------------------------------------------------------------------
# Small subprocess helper with error checking (no bare os.system).
# ---------------------------------------------------------------------------
def _run(cmd, cwd=None, env=None):
    log.info("RUN: %s", " ".join(str(c) for c in cmd))
    cp = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env,
                        capture_output=True, text=True)
    if cp.returncode != 0:
        log.error("Command failed (rc=%d): %s\nstdout:\n%s\nstderr:\n%s",
                  cp.returncode, " ".join(str(c) for c in cmd),
                  (cp.stdout or "")[-2000:], (cp.stderr or "")[-2000:])
        raise RuntimeError(f"Command failed (rc={cp.returncode}): {cmd[0]}")
    return cp


def preflight(require_openmm=True):
    """Verify AmberTools binaries (and optionally OpenMM) are importable/onPATH.
    Returns (ok: bool, problems: list[str])."""
    problems = []
    for t in REQUIRED_TOOLS:
        if shutil.which(t) is None:
            problems.append(f"Required tool not on PATH: {t}")
    if require_openmm:
        try:
            import openmm  # noqa: F401
        except Exception as e:  # noqa: BLE001
            problems.append(f"OpenMM not importable: {type(e).__name__}: {e}")
    return (not problems), problems


# ---------------------------------------------------------------------------
# Step 1: parametrize -> solvated complex prmtop/inpcrd (AmberTools/tleap).
# ---------------------------------------------------------------------------
def _ligand_to_pdb(ligand: Path, out_pdb: Path):
    """Convert an SDF/MOL2 ligand to PDB with explicit Hs (RDKit; obabel fallback)."""
    try:
        from rdkit import Chem
        if ligand.suffix.lower() == ".sdf":
            mol = next(iter(Chem.SDMolSupplier(str(ligand), removeHs=False)), None)
        else:
            mol = Chem.MolFromMol2File(str(ligand), removeHs=False)
        if mol is None:
            raise ValueError("RDKit could not parse ligand")
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.MolToPDBFile(mol, str(out_pdb))
        return out_pdb
    except Exception as e:  # noqa: BLE001
        log.warning("RDKit ligand->PDB failed (%s); trying obabel", e)
        if shutil.which("obabel") is None:
            raise RuntimeError("Cannot convert ligand to PDB: RDKit failed and "
                               "obabel not on PATH.") from e
        _run(["obabel", "-h", str(ligand), "-O", str(out_pdb)])
        return out_pdb


def split_ligand_file(ligand: Path, out_dir: Path) -> list[dict]:
    """Split a multi-record SDF/MOL2 into per-ligand single-record SDF files.

    Returns a list of {"index","name","sdf"} dicts (one per valid record).
    Names are sanitized for use as directory names. Reuses RDKit; mirrors the
    naming convention used by the Amber batch path (lig_<name>).
    """
    from rdkit import Chem
    out_dir.mkdir(parents=True, exist_ok=True)
    suppl = (Chem.SDMolSupplier(str(ligand), removeHs=False)
             if ligand.suffix.lower() == ".sdf"
             else _mol2_supplier(ligand))
    recs = []
    for i, mol in enumerate(suppl):
        if mol is None:
            log.warning("Ligand record %d failed to parse; skipping.", i)
            continue
        raw = (mol.GetProp("_Name") if (mol.HasProp("_Name")
               and mol.GetProp("_Name").strip()) else f"ligand_{i}")
        name = _sanitize(raw)
        sdf = out_dir / f"{name}.sdf"
        w = Chem.SDWriter(str(sdf))
        w.write(mol); w.close()
        recs.append({"index": i, "name": name, "sdf": sdf})
    if not recs:
        raise RuntimeError(f"No valid ligand records found in {ligand}")
    return recs


def _sanitize(name: str) -> str:
    import re as _re
    s = _re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return s.strip("_.") or "ligand"


def _mol2_supplier(path: Path):
    """Yield RDKit mols from a (possibly multi-record) MOL2 file."""
    from rdkit import Chem
    text = path.read_text(errors="ignore")
    marker = "@<TRIPOS>MOLECULE"
    if marker not in text:
        yield Chem.MolFromMol2File(str(path), removeHs=False)
        return
    blocks = [marker + b for b in text.split(marker)[1:]]
    for b in blocks:
        yield Chem.MolFromMol2Block(b, removeHs=False)



def _force_ligand_resname(pdb_path: Path, want: str,
                          protein_resnames: set | None = None) -> int:
    """Rewrite the ligand residue name in *pdb_path* to *want*.

    pdb4amber renames a GAFF ligand to UNL/<0>, which then will not match the
    antechamber-generated prep/frcmod (built with -rn <want>), causing tleap
    'Atom ... does not have a type' fatal errors. We rewrite any HETATM/ATOM
    residue that is clearly the small molecule (not a standard amino acid, not
    water/ion) to *want*. Returns the number of atom lines changed.
    """
    AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
          "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HID","HIE",
          "HIP","CYX","ASH","GLH","LYN","NME","ACE",
          # protonation/terminal variants
          "NALA","NARG","NASN","NASP","NCYS","NGLN","NGLU","NGLY","NHIS","NILE",
          "NLEU","NLYS","NMET","NPHE","NPRO","NSER","NTHR","NTRP","NTYR","NVAL",
          "CALA","CARG","CASN","CASP","CCYS","CGLN","CGLU","CGLY","CHIS","CILE",
          "CLEU","CLYS","CMET","CPHE","CPRO","CSER","CTHR","CTRP","CTYR","CVAL"}
    SOLV = {"HOH","WAT","TIP3","TP3","SPC","NA","NA+","CL","CL-","K","K+",
            "MG","CA","ZN","SO4","PO4"}
    keep = AA | SOLV
    if protein_resnames:
        keep = keep | {r.upper() for r in protein_resnames}
    lines = pdb_path.read_text().splitlines(keepends=True)
    out, n = [], 0
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 20:
            resn = ln[17:20].strip().upper()
            if resn not in keep:
                ln = ln[:17] + f"{want:>3}" + ln[20:]
                n += 1
        out.append(ln)
    pdb_path.write_text("".join(out))
    return n



def _pdb_atom_map_lines(prepi: Path) -> str:
    """Build tleap 'addPdbAtomMap' lines so case-different atom names in the
    loaded PDB (which PDB format upper-cases, e.g. 'CL1') match the antechamber
    prep names (which keep element case, e.g. 'Cl1').
    """
    try:
        lines = prepi.read_text().splitlines()
    except Exception:
        return ""
    names = []
    for ln in lines:
        toks = ln.split()
        if len(toks) >= 11 and toks[0].isdigit() and toks[1] not in ("DUMM",):
            nm = toks[1]
            if nm and nm.upper() != nm:
                names.append(nm)
    if not names:
        return ""
    pairs = "".join(f' {{ "{nm.upper()}" "{nm}" }}' for nm in dict.fromkeys(names))
    return f"addPdbAtomMap {{{pairs} }}\n"


def build_solvated_complex(protein: Path, ligand: Path, work: Path, *,
                           ligand_charge: int = 0, ligand_resname: str = "LIG",
                           charge_method: str = "bcc", buffer_A: float = 10.0,
                           env=None) -> dict:
    """Build solvated.prmtop/.inpcrd for the protein-ligand complex via tleap.

    Returns {'solvated_prmtop','solvated_inpcrd','complex_pdb'} paths.
    """
    work.mkdir(parents=True, exist_ok=True)

    # -- protein prep --
    prot_fixed = work / "protein_amber.pdb"
    _run(["pdb4amber", "-i", str(protein), "-o", str(prot_fixed)], cwd=work, env=env)

    # -- ligand prep: PDB -> mol2 (bcc charges) -> prepi + frcmod --
    lig_pdb = _ligand_to_pdb(Path(ligand), work / "ligand_h.pdb")
    lig_amber_pdb = work / "ligand_amber.pdb"
    _run(["pdb4amber", "-i", str(lig_pdb), "-o", str(lig_amber_pdb)], cwd=work, env=env)
    # pdb4amber may rename the ligand (e.g. to UNL/<0>); pin it to ligand_resname
    # NOW so both the antechamber prep (-rn) and the later tleap 'loadpdb lig'
    # refer to the same residue name (otherwise tleap: 'Atom ... has no type').
    _nlig = _force_ligand_resname(lig_amber_pdb, ligand_resname)
    log.info("Ligand-only PDB: set %d atom(s) to resname %s", _nlig, ligand_resname)
    lig_mol2 = work / "ligand.mol2"
    _run(["antechamber", "-fi", "pdb", "-i", str(lig_amber_pdb),
          "-fo", "mol2", "-o", str(lig_mol2),
          "-c", charge_method, "-nc", str(ligand_charge),
          "-rn", ligand_resname, "-pf", "y"], cwd=work, env=env)
    lig_prepi = work / "ligand.prepi"
    _run(["antechamber", "-i", str(lig_mol2), "-fi", "mol2",
          "-o", str(lig_prepi), "-fo", "prepi", "-pf", "y"], cwd=work, env=env)
    lig_frcmod = work / "ligand.frcmod"
    _run(["parmchk2", "-f", "prepi", "-i", str(lig_prepi),
          "-o", str(lig_frcmod)], cwd=work, env=env)

    # -- protein: strip ALL hydrogens so tleap rebuilds ff14SB-consistent H and
    #    assigns proper N-/C-terminal residues itself. The incoming protein
    #    (e.g. protein_omm.pdb) is typically already protonated by OpenMM with
    #    H names ff14SB cannot type (e.g. terminal MET carrying both 'H3' and
    #    'H', protonated GLU 'HE2'); -y removes them and tleap re-adds correct
    #    ones. This is the standard, robust Amber recipe. --
    prot_dry = work / "protein_dry.pdb"
    try:
        _run(["pdb4amber", "-y", "-i", str(prot_fixed), "-o", str(prot_dry)],
             cwd=work, env=env)
    except Exception:
        # Fallback: strip H lines ourselves (column 77-78 element == H).
        kept = []
        for ln in prot_fixed.read_text().splitlines(keepends=True):
            if ln.startswith(("ATOM", "HETATM")):
                elem = ln[76:78].strip().upper() if len(ln) >= 78 else ""
                aname = ln[12:16].strip()
                if elem == "H" or (not elem and aname[:1] == "H"):
                    continue
            kept.append(ln)
        prot_dry.write_text("".join(kept))

    # Keep complex.pdb (heavy-atom protein + parameterized ligand) for reference
    # and for MMPBSA topology generation downstream.
    complex_amber_pdb = work / "complex_amber.pdb"
    complex_amber_pdb.write_text(prot_dry.read_text() + lig_amber_pdb.read_text())

    # -- tleap: load protein and ligand as SEPARATE units, combine, then solvate.
    #    Loading the ligand from its own untouched PDB keeps atom names (e.g.
    #    'CL1') matching the antechamber prep, and loading the protein heavy
    #    atoms lets ff14SB add correctly-named H + termini. --
    solv_prm = work / "solvated.prmtop"
    solv_crd = work / "solvated.inpcrd"
    tleap_in = work / "complex.tleap.in"
    atom_map = _pdb_atom_map_lines(lig_prepi)   # fix Cl1/CL1 element-case mismatch
    if atom_map:
        log.info("tleap atom-name map for ligand: %s", atom_map.strip())
    tleap_in.write_text(
        "source leaprc.protein.ff14SB\n"
        "source leaprc.gaff\n"
        "source leaprc.water.tip3p\n"
        f"loadamberprep {lig_prepi.name}\n"
        f"loadamberparams {lig_frcmod.name}\n"
        + atom_map +
        f"prot = loadpdb {prot_dry.name}\n"
        f"lig = loadpdb {lig_amber_pdb.name}\n"
        "mol = combine {prot lig}\n"
        f"solvatebox mol TIP3PBOX {buffer_A}\n"
        "addions mol Na+ 0\n"
        "addions mol Cl- 0\n"
        f"saveamberparm mol {solv_prm.name} {solv_crd.name}\n"
        f"savepdb mol {complex_amber_pdb.name}\n"
        "quit\n")
    tleap_out = work / "complex.tleap.out"
    cp = subprocess.run(["tleap", "-s", "-f", tleap_in.name],
                        cwd=str(work), env=env, capture_output=True, text=True)
    tleap_out.write_text((cp.stdout or "") + "\n----STDERR----\n" + (cp.stderr or ""))
    if not (solv_prm.exists() and solv_crd.exists()):
        tail = ""
        try:
            _all = tleap_out.read_text().splitlines()
            fatal = [l for l in _all if "FATAL" in l or "Fatal" in l]
            tail = "\n".join((fatal[:8] or _all[-25:]))
        except Exception:
            pass
        log.error("tleap did not produce topology. See %s\n--- tleap errors ---\n%s",
                  tleap_out, tail)
        raise RuntimeError(
            f"tleap failed to build solvated complex (see {tleap_out}).\n{tail}")
    log.info("Built solvated complex: %s", solv_prm)
    return {"solvated_prmtop": solv_prm, "solvated_inpcrd": solv_crd,
            "complex_pdb": complex_amber_pdb}


# ---------------------------------------------------------------------------
# Step 2: OpenMM MD -> NetCDF trajectory matching solvated.prmtop.
# ---------------------------------------------------------------------------
def run_openmm_md(solvated_prmtop: Path, solvated_inpcrd: Path, work: Path, *,
                  prod_ns: float = 10.0, equil_ns: float = 1.0,  # final51: unified 10 ns (both engines)
                  timestep_ps: float = 0.004, temperature_K: float = 300.0,
                  platform_name: str = "CUDA", report_ps: float = 10.0) -> Path:
    """Run minimize/equilibrate/produce in OpenMM. Returns the NetCDF traj path.

    Frames are written for the PRODUCTION phase only. dt=4 fs uses HMR
    (hydrogenMass=1.5 amu) + HBond constraints, matching the reference recipe.
    """
    import openmm
    from openmm import app, unit

    prmtop = app.AmberPrmtopFile(str(solvated_prmtop))
    inpcrd = app.AmberInpcrdFile(str(solvated_inpcrd))

    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        hydrogenMass=1.5 * unit.amu,
        ewaldErrorTolerance=0.0005)
    system.addForce(openmm.MonteCarloBarostat(
        1.0 * unit.atmospheres, temperature_K * unit.kelvin, 25))

    dt = timestep_ps * unit.picoseconds
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin, 1.0 / unit.picosecond, dt)

    platform = openmm.Platform.getPlatformByName(platform_name)
    sim = app.Simulation(prmtop.topology, system, integrator, platform)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    steps_per_ns = int(round(1000.0 / timestep_ps))   # ps per ns / dt(ps)
    equil_steps = int(equil_ns * steps_per_ns)
    prod_steps = int(prod_ns * steps_per_ns)
    report_interval = max(1, int(report_ps / timestep_ps))

    log.info("OpenMM MD: equil=%d steps, prod=%d steps, dt=%.3f ps, platform=%s",
             equil_steps, prod_steps, timestep_ps, platform_name)

    log.info("Minimizing...")
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    if equil_steps > 0:
        log.info("Equilibrating %d steps...", equil_steps)
        sim.step(equil_steps)

    # NetCDF trajectory for the production phase (Amber-native; MMPBSA.py -y).
    traj = work / "production.nc"
    try:
        from openmm.app import NetCDFReporter  # OpenMM >=7.7 (parmed-backed)
        sim.reporters.append(NetCDFReporter(str(traj), report_interval))
    except Exception:  # pragma: no cover - fall back to mdtraj reporter
        try:
            from mdtraj.reporters import NetCDFReporter as MDTrajNetCDF
            sim.reporters.append(MDTrajNetCDF(str(traj), report_interval))
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(
                "No NetCDF reporter available (need OpenMM>=7.7 with NetCDF "
                "support or mdtraj installed).") from e

    sim.reporters.append(app.StateDataReporter(
        str(work / "md_log.txt"), report_interval, step=True,
        potentialEnergy=True, temperature=True, speed=True,
        totalSteps=prod_steps, separator="\t"))

    log.info("Producing %d steps -> %s", prod_steps, traj)
    sim.step(prod_steps)
    if not traj.exists():
        raise RuntimeError(f"OpenMM finished but no trajectory at {traj}")
    return traj


# ---------------------------------------------------------------------------
# Step 3: score (reuse existing TopologySplitter + MMGBSAAnalyzer).
# ---------------------------------------------------------------------------
def score_mmgbsa(solvated_prmtop: Path, trajectory: Path, work: Path, *,
                 igb: int = 8, salt: float = 0.15, ligand_resname: str = "LIG",
                 interval: int = 2, decompose: bool = False) -> Path:
    from .topology import TopologySplitter
    from .mmgbsa import MMGBSAAnalyzer

    # GB radius set MUST match igb or mmpbsa_py_energy fails on complex.prmtop:
    #   igb 1 -> mbondi, igb 2/5 -> mbondi2, igb 7 -> bondi, igb 8 -> mbondi3.
    _IGB_RADII = {1: "mbondi", 2: "mbondi2", 5: "mbondi2", 7: "bondi", 8: "mbondi3"}
    radii = _IGB_RADII.get(int(igb), "mbondi3")
    splitter = TopologySplitter(work, ligand_resname=ligand_resname, radii=radii)
    # force=True so a previously-built topo/ with the WRONG radii is rebuilt.
    topos = splitter.split(solvated_prmtop, force=True)
    TopologySplitter.sanity_check(topos)
    topo_dir = topos["complex"].parent  # all 4 prmtops live here (work/topo)

    analyzer = MMGBSAAnalyzer(
        topo_dir, trajectory, work / "mmgbsa",
        igb=igb, saltcon=salt, interval=interval, decompose=decompose)
    return analyzer.run()


# ---------------------------------------------------------------------------
# Orchestration + CLI.
# ---------------------------------------------------------------------------
def run_batch(args) -> int:
    """Multi-ligand mode: split the ligand file, run one MM-GBSA per ligand into
    <workdir>/lig_<name>/, then (best-effort) aggregate+rank via batch_aggregate.

    Layout matches amber_md.batch_aggregate's contract:
      <workdir>/lig_<name>/mmgbsa/FINAL_RESULTS_MMPBSA.dat
    """
    work = Path(args.workdir).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)

    ok, problems = preflight(require_openmm=not args.dry_run)
    for p in problems:
        log.error("PREFLIGHT: %s", p)
    if not ok and not args.dry_run:
        log.error("Preflight failed; aborting batch.")
        return 2

    recs = split_ligand_file(Path(args.ligand), work / "_ligands")
    log.info("Batch MM-GBSA: %d ligand(s) -> %s", len(recs), work)

    failures = 0
    for r in recs:
        lig_dir = work / f"lig_{r['name']}"
        lig_dir.mkdir(parents=True, exist_ok=True)
        log.info("=== Ligand %s (record %d) -> %s ===",
                 r["name"], r["index"], lig_dir)
        if args.dry_run:
            log.info("--dry-run: would run MM-GBSA for %s", r["name"])
            continue
        try:
            # Delegate to run_pipeline so the per-ligand run gets the SAME
            # behaviour as a single run -- including the resume guard (skip
            # build+MD when production.nc already exists) and every analysis
            # fix. We clone args and point it at this ligand's dir + sdf.
            import copy as _copy
            la = _copy.copy(args)
            la.batch = False
            la.workdir = str(lig_dir)
            la.ligand = str(r["sdf"])
            rc = run_pipeline(la)
            if rc != 0:
                failures += 1
                log.error("Ligand %s: run_pipeline returned rc=%d", r["name"], rc)
            else:
                log.info("Ligand %s: DONE", r["name"])
        except Exception as e:  # noqa: BLE001
            failures += 1
            log.exception("Ligand %s FAILED: %s", r["name"], e)

    # Best-effort aggregate + rank (reuses the Amber batch ranking).
    if not args.dry_run:
        try:
            import subprocess, sys as _sys
            subprocess.run([_sys.executable, "-m", "amber_md.batch_aggregate",
                            str(work)], check=False)
            log.info("Aggregation attempted; see %s/INDEX.html", work)
        except Exception as e:  # noqa: BLE001
            log.warning("batch_aggregate failed (%s); per-ligand .dat files "
                        "are still in lig_*/mmgbsa/", e)

    log.info("Batch complete: %d ligand(s), %d failure(s).",
             len(recs), failures)
    return 1 if failures == len(recs) else 0


def run_pipeline(args) -> int:
    work = Path(args.workdir).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)

    ok, problems = preflight(require_openmm=not args.dry_run)
    for p in problems:
        log.error("PREFLIGHT: %s", p)
    if not ok and not args.dry_run:
        log.error("Preflight failed; aborting. (Install AmberTools + OpenMM in "
                  "the env, or run with --dry-run to validate inputs only.)")
        return 2

    plan = {
        "protein": str(Path(args.protein).resolve()),
        "ligand": str(Path(args.ligand).resolve()),
        "workdir": str(work),
        "prod_ns": args.prod_ns, "equil_ns": args.equil_ns,
        "igb": args.igb, "salt": args.salt, "platform": args.platform,
        "ligand_charge": args.ligand_charge, "ligand_resname": args.ligand_resname,
    }
    log.info("OpenMM MM-GBSA plan: %s", plan)
    if args.dry_run:
        log.info("--dry-run: inputs validated, plan printed; not executing.")
        return 0

    # ---- resume guard (default ON) -----------------------------------------
    # Building the topology + running 5 ns of MD is the expensive part (~15 min
    # on a GPU). If a previous run already produced the solvated topology and a
    # production trajectory, skip straight to scoring. This makes re-runs after
    # an analysis-stage failure (or a deliberate re-score) near-instant, and
    # lets a crashed batch member resume instead of repeating MD.
    # Pass --no-resume to force a clean rebuild + fresh MD.
    solv_prm = work / "solvated.prmtop"
    solv_crd = work / "solvated.inpcrd"
    traj = work / "production.nc"
    resume = getattr(args, "resume", True)

    def _nonempty(p):
        return p.exists() and p.stat().st_size > 0

    have_build = _nonempty(solv_prm) and _nonempty(solv_crd)
    have_traj = _nonempty(traj)

    if resume and have_build and have_traj:
        log.info("RESUME: found existing topology (%s) and trajectory (%s); "
                 "skipping build + MD, scoring only. Use --no-resume to "
                 "force a fresh run.", solv_prm.name, traj.name)
        built = {"solvated_prmtop": solv_prm, "solvated_inpcrd": solv_crd}
    else:
        if resume and (have_build or have_traj):
            log.info("RESUME: partial artifacts present (build=%s, traj=%s); "
                     "rebuilding the missing stage(s).", have_build, have_traj)
        built = build_solvated_complex(
            Path(args.protein), Path(args.ligand), work,
            ligand_charge=args.ligand_charge, ligand_resname=args.ligand_resname,
            charge_method=args.charge_method, buffer_A=args.buffer)
        if resume and have_traj and _nonempty(built["solvated_prmtop"]):
            log.info("RESUME: reusing existing trajectory %s.", traj.name)
        else:
            traj = run_openmm_md(
                built["solvated_prmtop"], built["solvated_inpcrd"], work,
                prod_ns=args.prod_ns, equil_ns=args.equil_ns,
                temperature_K=args.temperature, platform_name=args.platform)

    dat = score_mmgbsa(
        built["solvated_prmtop"], traj, work,
        igb=args.igb, salt=args.salt, ligand_resname=args.ligand_resname,
        interval=args.interval, decompose=args.decomp)
    # final63: drop a self-describing engine marker next to the result so
    # Results-Compare / results_lib can label the row "OpenMM" unambiguously
    # (the Amber pmemd path never writes mmgbsa/engine.json).
    try:
        mm_dir = work / "mmgbsa"
        mm_dir.mkdir(parents=True, exist_ok=True)
        (mm_dir / "engine.json").write_text(json.dumps({
            "engine": "OpenMM", "method": "MM-GBSA",
            "runner": "amber_md.mmgbsa_openmm", **plan}, indent=2))
    except Exception as _e:  # noqa: BLE001
        log.warning("could not write engine marker: %s", _e)
    log.info("DONE. MM-GBSA result: %s", dat)
    return 0



def _submit_to_lsf(a, argv) -> int:
    """Write a self-contained LSF (#BSUB) script that re-runs THIS MM-GBSA job
    on the GPU queue with --submit local, then `bsub <` it.

    Reuses the cluster conventions from amber_md.submit / config.HPCConfig
    (queue, project, walltime, module loads, env activation). Returns 0 on a
    successful submission, non-zero otherwise. The science code is unchanged --
    the compute node simply runs the same CLI locally.
    """
    import shlex as _shlex
    import subprocess as _sp
    work = Path(a.workdir).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)

    # Defaults from HPCConfig where the user did not override.
    try:
        from .config import HPCConfig
        hpc = HPCConfig()
    except Exception:
        hpc = None
    queue = a.queue or (getattr(hpc, "queue_gpu", None) or "gpu")
    project = a.project or (getattr(hpc, "project", None) or "")
    if not project:
        log.warning("No LSF project (-P) resolved; omitting it. The GPU "
                    "queue may reject the job -- pass --project "
                    "your-project or set HPCConfig.project.")
    walltime = a.walltime or (getattr(hpc, "walltime", None) or "24:00")
    # IMPORTANT: the OpenMM MM-GBSA pipeline is fully self-contained inside the
    # conda env (its own tleap/antechamber/MMPBSA.py/parmed). Loading the system
    # `amber/*` module prepends /share/apps/amber/.../bin to PATH and shadows
    # those env tools with an OLD parmed that does `from numpy.compat import ...`
    # -- removed in NumPy >= 1.24 -- crashing ante-MMPBSA.py / parmed with
    # "No module named 'numpy.compat'". So we DROP any amber* module here and
    # rely on the env. CUDA/gcc are still loaded for the GPU runtime.
    _all_modules = list(getattr(hpc, "modules", ()) or [])
    if a.modules is not None:
        _all_modules = [x for x in a.modules.split(",") if x.strip()]
    modules = [x for x in _all_modules if not x.lower().startswith("amber")]
    _dropped = [x for x in _all_modules if x.lower().startswith("amber")]
    n_gpu = int(a.n_gpu or getattr(hpc, "n_gpu", 1) or 1)
    n_cpu = int(getattr(hpc, "n_cpu", 8) or 8)

    # Rebuild the original argv but force --submit local so the node runs it.
    # Also make the resume choice EXPLICIT in the inner command. This is what
    # makes a requeued / crashed-then-resubmitted job pick up where it left off:
    # the node re-runs in the SAME workdir (cd below) and, with --resume, skips
    # the expensive build + 5 ns MD when solvated.prmtop + production.nc already
    # exist. We strip any pre-existing --submit/--resume/--no-resume tokens and
    # re-append a single canonical pair so the behaviour is deterministic even
    # if argparse defaults change later.
    inner = list(argv) if argv is not None else sys.argv[1:]
    out = []
    skip = False
    for i, tok in enumerate(inner):
        if skip:
            skip = False
            continue
        if tok == "--submit":
            skip = True            # drop "--submit gpu"
            continue
        if tok.startswith("--submit=") or tok in ("--resume", "--no-resume"):
            continue               # drop; we re-add the canonical resume flag
        out.append(tok)
    out += ["--submit", "local"]
    # Honour the user's actual resume choice (defaults True). On the GPU/requeue
    # path resume is what enables crash recovery, so it is ON unless the user
    # explicitly passed --no-resume.
    resume = getattr(a, "resume", True)
    out += ["--resume"] if resume else ["--no-resume"]
    py = sys.executable
    inner_cmd = " ".join(_shlex.quote(x) for x in [py, "-m", "amber_md.mmgbsa_openmm", *out])

    # Conda activation: prefer an explicit env, else the running interpreter's.
    conda_env = a.conda_env
    if not conda_env:
        # sys.prefix points at the active env (…/envs/<name>)
        conda_env = sys.prefix
    activate = (f'source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" '
                f'2>/dev/null || true\n'
                f'conda activate {_shlex.quote(conda_env)} 2>/dev/null || '
                f'source activate {_shlex.quote(conda_env)} 2>/dev/null || true\n')

    # Defensive: remove any system-amber bin dirs from PATH so the env's
    # AmberTools always win, even if a site profile re-injected them.
    path_clean = (
        'export PATH="$(printf "%s" "$PATH" | tr ":" "\n" '
        '| grep -v -E "/share/apps/amber|/amber[0-9]*/bin" | paste -sd: -)"\n'
    )
    if _dropped:
        log.info("Dropping system module(s) %s for OpenMM MM-GBSA (env provides "
                 "AmberTools; avoids parmed/numpy.compat crash).", _dropped)

    workflow_root = str(Path(__file__).resolve().parent.parent)
    job_name = f"mmgbsa_{work.name}"
    log_prefix = str(work / "lsf_mmgbsa")

    script = "#!/bin/bash\n"
    script += f"#BSUB -q {queue}\n"
    if project:
        script += f"#BSUB -P {project}\n"
    script += f"#BSUB -J {job_name}\n#BSUB -W {walltime}\n"
    script += f"#BSUB -o {log_prefix}.%J.out\n#BSUB -e {log_prefix}.%J.err\n"
    script += f"#BSUB -n {n_cpu}\n"
    script += '#BSUB -R "span[hosts=1]"\n'
    script += '#BSUB -R "rusage[mem=4096]"\n'
    script += f'#BSUB -gpu "num={n_gpu}"\n'
    # Auto-requeue on transient/abnormal termination (node failure, preemption,
    # etc.). When LSF requeues the job it re-runs THIS script in the SAME
    # workdir; with --resume the node then skips completed build+MD and resumes
    # from production.nc. Exit code 0 (success) is NOT requeued. Disable with
    # --no-requeue. We request requeue for the common transient exit codes plus
    # the LSF preemption/limit signals via -r as a fallback.
    requeue = getattr(a, "requeue", True)
    max_requeue = int(getattr(a, "max_requeue", 3) or 3)
    if requeue:
        # -Q "all ~0" => requeue on any non-zero exit except 0; safe because a
        # finished run returns 0 and a resumed run is cheap (scores only).
        # We bound the total attempts (below) so a DETERMINISTIC failure (a real
        # bug that fails identically every time) cannot loop forever.
        script += '#BSUB -Q "all ~0"\n'
        script += "#BSUB -r\n"            # rerunnable: requeue on host failure
    script += "\n"
    script += "module purge 2>/dev/null || true\n"
    for mod in modules:
        script += f"module load {mod} 2>/dev/null || true\n"
    script += activate
    script += path_clean
    script += f'export PYTHONPATH="{workflow_root}${{PYTHONPATH:+:$PYTHONPATH}}"\n'
    script += "set -euo pipefail\n\n"
    script += f"cd {_shlex.quote(str(work))}\n"
    script += f"echo \"[mmgbsa] node=$(hostname) gpu=${{CUDA_VISIBLE_DEVICES:-?}}\"\n"
    if requeue:
        # Bounded auto-requeue: track attempts in the workdir. If we exceed the
        # cap, exit with code 0 in a way that does NOT trigger another requeue
        # (the -Q set excludes 0) and leave a marker so the failure is visible.
        # This prevents an infinite resubmit loop on a deterministic failure
        # while still allowing genuine transient retries (node death, preempt).
        script += (
            f"ATT_FILE={_shlex.quote(str(work / '.requeue_attempts'))}\n"
            f"ATT=$(cat \"$ATT_FILE\" 2>/dev/null || echo 0)\n"
            f"ATT=$((ATT+1)); echo \"$ATT\" > \"$ATT_FILE\"\n"
            f"echo \"[mmgbsa] requeue attempt $ATT of {max_requeue}\"\n"
            f"if [ \"$ATT\" -gt {max_requeue} ]; then\n"
            f"  echo \"[mmgbsa] exceeded max_requeue={max_requeue}; not retrying.\" >&2\n"
            f"  echo \"requeue_cap_exceeded after $ATT attempts\" "
            f"> {_shlex.quote(str(work / 'REQUEUE_GIVING_UP.txt'))}\n"
            f"  exit 0\n"
            f"fi\n"
        )
    script += inner_cmd + "\n"
    if requeue:
        # On SUCCESS, clear the attempt counter so a later deliberate re-submit
        # of the same workdir starts fresh.
        script += (
            f"rm -f {_shlex.quote(str(work / '.requeue_attempts'))} "
            f"{_shlex.quote(str(work / 'REQUEUE_GIVING_UP.txt'))} 2>/dev/null || true\n"
        )

    script_path = work / "submit_mmgbsa.lsf"
    script_path.write_text(script)
    log.info("Wrote LSF script: %s", script_path)

    if shutil.which("bsub") is None:
        log.error("bsub not found on PATH -- cannot submit. The script is ready "
                  "at %s; submit it manually with: bsub < %s",
                  script_path, script_path)
        return 3
    try:
        cp = _sp.run(["bsub"], stdin=open(script_path), capture_output=True,
                     text=True, check=False)
    except Exception as e:  # noqa: BLE001
        log.error("bsub invocation failed: %s", e)
        return 3
    sys.stdout.write(cp.stdout)
    sys.stderr.write(cp.stderr)
    if cp.returncode != 0:
        log.error("bsub returned %d", cp.returncode)
        return cp.returncode
    log.info("Submitted MM-GBSA to LSF queue '%s' (project '%s', walltime %s).",
             queue, project, walltime)
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="python -m amber_md.mmgbsa_openmm",
        description="OpenMM-MD + AmberTools MMPBSA.py MM-GBSA (Interpretation A).")
    p.add_argument("--protein", required=True, help="Protein PDB")
    p.add_argument("--ligand", required=True, help="Ligand SDF/MOL2 (single; or multi-record with --batch)")
    p.add_argument("--workdir", required=True)
    p.add_argument("--prod-ns", type=float, default=10.0, dest="prod_ns")  # final51: unified 10 ns
    p.add_argument("--equil-ns", type=float, default=1.0, dest="equil_ns")
    p.add_argument("--temperature", type=float, default=300.0)
    p.add_argument("--platform", default="CUDA", choices=["CUDA", "OpenCL", "CPU"])
    p.add_argument("--igb", type=int, default=8, choices=[1, 2, 5, 7, 8])
    p.add_argument("--salt", type=float, default=0.15)
    p.add_argument("--interval", type=int, default=2,
                   help="MMPBSA frame stride")
    p.add_argument("--ligand-charge", type=int, default=0, dest="ligand_charge")
    p.add_argument("--ligand-resname", default="LIG", dest="ligand_resname")
    p.add_argument("--charge-method", default="bcc",
                   choices=["bcc", "gas"], dest="charge_method")
    p.add_argument("--buffer", type=float, default=10.0,
                   help="Solvent box buffer (Angstrom)")
    p.add_argument("--decomp", action="store_true",
                   help="Per-residue decomposition")
    p.add_argument("--batch", action="store_true",
                   help="Treat --ligand as multi-record; one job per ligand.")
    p.add_argument("--dry-run", action="store_true", dest="dry_run")
    # ---- execution target -------------------------------------------------
    p.add_argument("--submit", default="local", choices=["local", "gpu"],
                   help="local = run here (default); gpu = submit to the LSF "
                        "GPU queue via bsub and exit.")
    p.add_argument("--queue", default="gpu",
                   help="LSF GPU queue name when --submit gpu (default: gpu).")
    p.add_argument("--walltime", default="24:00",
                   help="LSF walltime (HH:MM) when --submit gpu.")
    p.add_argument("--project", default=None,
                   help="LSF project (-P). Defaults to HPCConfig.project.")
    p.add_argument("--conda-env", default=None, dest="conda_env",
                   help="Conda env to activate inside the LSF job (name or "
                        "path). Defaults to the env of the current python.")
    p.add_argument("--n-gpu", type=int, default=1, dest="n_gpu",
                   help="GPUs to request when --submit gpu (default: 1).")
    p.add_argument("--resume", dest="resume", action="store_true", default=True,
                   help="Reuse an existing solvated topology + production.nc in "
                        "the workdir and skip build+MD, scoring only (DEFAULT). "
                        "Makes re-runs after an analysis failure near-instant.")
    p.add_argument("--no-resume", dest="resume", action="store_false",
                   help="Force a clean rebuild and a fresh MD run, ignoring any "
                        "existing artifacts in the workdir.")
    p.add_argument("--requeue", dest="requeue", action="store_true", default=True,
                   help="On --submit gpu, add LSF auto-requeue (#BSUB -Q/-r) so "
                        "a crashed/preempted job is resubmitted and resumes from "
                        "production.nc in the same workdir (DEFAULT).")
    p.add_argument("--no-requeue", dest="requeue", action="store_false",
                   help="Disable LSF auto-requeue for the GPU job.")
    p.add_argument("--max-requeue", dest="max_requeue", type=int, default=3,
                   help="Max LSF auto-requeue attempts before giving up, to "
                        "avoid an infinite resubmit loop on a deterministic "
                        "failure (default: 3).")
    p.add_argument("--modules", default=None,
                   help="Comma-separated `module load` list for the LSF job. "
                        "Overrides HPCConfig.modules. Any 'amber*' module is "
                        "dropped automatically (the conda env provides its own "
                        "AmberTools; the system module breaks parmed/numpy).")
    a = p.parse_args(argv)

    if a.submit == "gpu":
        return _submit_to_lsf(a, argv)
    return run_batch(a) if a.batch else run_pipeline(a)


if __name__ == "__main__":
    sys.exit(main() or 0)
