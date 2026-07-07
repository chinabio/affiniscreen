"""Multi-ligand batch screening (v2.4.7 - threads --decomp through to each ligand)."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import os, sys, re, time, json, subprocess, argparse
from pathlib import Path
from datetime import datetime
from .logger import get_logger
log = get_logger()

def _have_rdkit():
    try:
        import rdkit  # noqa
        return True
    except ImportError:
        return False

def _extract_charge_from_sdf(text):
    for tag in ("i_epik_Tot_Q", "i_user_NetCharge", "r_user_NetCharge",
                "TOTAL_CHARGE", "CHARGE", "NetCharge"):
        m = re.search(rf"^>\s*<{re.escape(tag)}>[^\n]*\n\s*(-?\d+(?:\.\d+)?)",
                      text, re.MULTILINE)
        if m:
            try: return int(float(m.group(1).strip()))
            except Exception: pass
    return None

def _is_v3000(text):
    head = "\n".join(text.splitlines()[:10])
    return "V3000" in head and "V2000" not in head

def preflight_protein(protein_pdb, scratch):
    """Run tleap on the protein alone. Returns (ok: bool, message: str)."""
    scratch = Path(scratch); scratch.mkdir(parents=True, exist_ok=True)
    test_in = scratch / "_preflight_protein.in"
    test_in.write_text(
        "source leaprc.protein.ff19SB\n"
        f"p = loadpdb {protein_pdb}\n"
        "charge p\n"
        "quit\n"
    )
    try:
        cp = subprocess.run(["tleap", "-f", str(test_in)],
                            capture_output=True, text=True, timeout=120,
                            cwd=str(scratch))
    except Exception as e:
        return False, f"tleap invocation failed: {e}"
    out = cp.stdout + cp.stderr
    fatal = [ln for ln in out.splitlines() if "FATAL" in ln]
    if fatal:
        msg = (f"tleap rejected {protein_pdb}:\n"
               + "\n".join("  " + l for l in fatal[:20])
               + (f"\n  ...({len(fatal)-20} more FATAL lines)" if len(fatal)>20 else ""))
        return False, msg
    m = re.search(r"Total unperturbed charge:\s*(-?\d+\.\d+)", out)
    chg = m.group(1) if m else "?"
    return True, f"tleap OK on {Path(protein_pdb).name}, net charge = {chg}"

def _split_with_rdkit(src, out_dir):
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    log.info("Splitting %s via RDKit (handles V3000, extracts charges)...", src.name)
    suppl = Chem.SDMolSupplier(str(src), removeHs=False, sanitize=False)
    out = []
    for i, mol in enumerate(suppl, 1):
        if mol is None:
            log.warning("  record #%d: RDKit could not parse -- skipping", i); continue
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception as e:
            log.warning("  record #%d: sanitize warning: %s", i, e)
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else f"MOL{i}"
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", name)[:40] or f"MOL{i}"
        charge = None
        for tag in ("i_epik_Tot_Q", "i_user_NetCharge", "TOTAL_CHARGE", "CHARGE"):
            if mol.HasProp(tag):
                try: charge = int(float(mol.GetProp(tag))); break
                except Exception: pass
        if charge is None:
            charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        fname = out_dir / f"lig_{i:04d}_{safe}.sdf"
        w = Chem.SDWriter(str(fname)); w.SetForceV3000(False); w.write(mol); w.close()
        out.append((i, safe, fname, charge))
    return out

def _split_multi_mol2(src, out_dir):
    text = src.read_text(errors="replace").replace("\r\n","\n").replace("\r","\n")
    parts = text.split("@<TRIPOS>MOLECULE")
    out = []
    for i, body in enumerate(parts[1:], 1):
        lines = body.splitlines()
        name = "MOL"
        for ln in lines:
            s = ln.strip()
            if s:
                name = re.sub(r"[^A-Za-z0-9_]+", "_", s)[:40] or "MOL"; break
        fname = out_dir / f"lig_{i:04d}_{name}.mol2"
        fname.write_text("@<TRIPOS>MOLECULE" + body)
        out.append((i, name, fname, None))
    return out

def _split_multi_sdf_legacy(src, out_dir):
    text = src.read_text(errors="replace").replace("\r\n","\n").replace("\r","\n")
    records = re.split(r"^\$\$\$\$\s*$\n?", text, flags=re.MULTILINE)
    out = []; counter = 0
    for body in records:
        body = body.rstrip("\n").rstrip()
        if not body.strip(): continue
        if _is_v3000(body):
            log.error("Legacy splitter cannot handle V3000."); return []
        counter += 1
        lines = body.splitlines()
        name = "MOL"
        for ln in lines:
            s = ln.strip()
            if s:
                name = re.sub(r"[^A-Za-z0-9_]+", "_", s)[:40] or "MOL"; break
        charge = _extract_charge_from_sdf(body)
        fname = out_dir / f"lig_{counter:04d}_{name}.sdf"
        fname.write_text(body + "\n$$$$\n")
        out.append((counter, name, fname, charge))
    return out

def _discover_dir(src):
    out = []; counter = 0
    for f in sorted(src.iterdir()):
        if not f.is_file(): continue
        ext = f.suffix.lower()
        if ext not in {".mol2",".sdf",".pdb",".mol"}: continue
        counter += 1
        name = re.sub(r"[^A-Za-z0-9_]+","_", f.stem)[:40] or f"MOL{counter}"
        charge = None
        if ext == ".sdf":
            try:
                text = f.read_text(errors="replace")
                if _is_v3000(text):
                    log.error("File %s is V3000 -- skipping. Use prep_ligands.py.", f.name); continue
                charge = _extract_charge_from_sdf(text)
            except Exception as e:
                log.warning("Could not parse charge from %s: %s", f.name, e)
        out.append((counter, name, f, charge))
    return out

def discover_ligands(input_path, scratch):
    p = Path(input_path).resolve()
    if p.is_dir():
        log.info("Ligands: scanning directory %s", p)
        return _discover_dir(p)
    if not p.exists():
        raise FileNotFoundError(f"Ligand input not found: {p}")
    ext = p.suffix.lower()
    if ext == ".sdf":
        text = p.read_text(errors="replace")
        if _is_v3000(text) and not _have_rdkit():
            log.error("="*70)
            log.error("Input is V3000 but RDKit not installed.")
            log.error("Run: python prep_ligands.py %s ./cleaned/", p.name)
            log.error("="*70)
            raise SystemExit(1)
        split_dir = scratch / "split_sdf"; split_dir.mkdir(parents=True, exist_ok=True)
        mols = (_split_with_rdkit(p, split_dir) if _have_rdkit()
                else _split_multi_sdf_legacy(p, split_dir))
        log.info("Split %s into %d SDF record(s)", p.name, len(mols))
        return mols
    if ext == ".mol2":
        split_dir = scratch / "split_mol2"; split_dir.mkdir(parents=True, exist_ok=True)
        mols = _split_multi_mol2(p, split_dir)
        log.info("Split %s into %d MOL2 record(s)", p.name, len(mols))
        return mols
    log.info("Treating %s as a single ligand", p)
    return [(1, p.stem, p, None)]

def _count_active_jobs(user, queue=None):
    cmd = ["bjobs","-u",user,"-noheader"]
    if queue: cmd += ["-q", queue]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        if cp.returncode != 0: return 0
        return sum(1 for l in cp.stdout.splitlines()
                   if l.strip() and "No unfinished" not in l)
    except Exception as e:
        log.warning("Could not count active jobs: %s", e); return 0

def _wait_for_capacity(user, queue, cap, poll=30):
    while True:
        n = _count_active_jobs(user, queue)
        if n < cap:
            log.info("LSF capacity OK: %d/%d in queue.", n, cap); return
        log.info("LSF at capacity: %d/%d active; sleeping %ds...", n, cap, poll)
        time.sleep(poll)

def _submit_one(idx, name, lig_file, charge, args, batch_dir):
    workdir = batch_dir / f"lig_{idx:04d}_{name}"
    workdir.mkdir(parents=True, exist_ok=True)
    log_file = workdir / "pipeline.log"
    effective_charge = charge if charge is not None else args.ligand_charge
    cmd = [sys.executable, str(args.run_amber_py),
           "--protein-file", str(args.protein_file),
           "--ligand-file",  str(lig_file),
           "--lig-resname",  args.lig_resname,
           "--ligand-charge", str(effective_charge),
           "--charge-method", args.charge_method,
           "--workdir", str(workdir),
           "--prod-ns", str(args.prod_ns),
           "--equil-ns", str(args.equil_ns),
           "--walltime", args.walltime,
           "--job-name", f"{args.job_prefix}_{idx:04d}"[:24],
           "--n-gpu", str(args.n_gpu),
           "--salt", str(args.salt),
           "--ion-method", args.ion_method,
           "--no-monitor"]
    if args.project: cmd += ["--project", args.project]
    if args.queue:   cmd += ["--queue", args.queue]
    if args.no_gbsa: cmd.append("--no-gbsa")
    if not args.protonation: cmd.append("--no-protonation")
    # v2.4.7: forward decomposition flags
    if args.decomp:
        cmd.append("--decomp")
        if args.decomp_residues:
            cmd += ["--decomp-residues", args.decomp_residues]
    log.info("Launching #%04d %-30s (charge=%+d) -> %s",
             idx, name, effective_charge, workdir.name)
    with open(log_file, "w") as f:
        f.write(f"# Ligand #{idx}: {name}  charge={effective_charge}\n")
        f.write(f"# Source: {lig_file}\n# Launched: {datetime.now().isoformat()}\n")
        f.write(f"# Command: {' '.join(cmd)}\n\n"); f.flush()
        proc = subprocess.Popen(cmd, cwd=str(args.run_amber_py.parent),
                                stdout=f, stderr=subprocess.STDOUT,
                                start_new_session=True)
    return proc, workdir, log_file

def _parse_jobid(log_file, wait):
    """v2.4.3: requires BOTH 'Pipeline finished' AND 'Job <NNN>' for success."""
    deadline = time.time() + wait
    while time.time() < deadline:
        if log_file.exists():
            text = log_file.read_text()
            has_jobid   = re.search(r"Job <(\d+)>", text)
            has_done    = "Pipeline finished" in text
            has_crash   = ("Traceback (most recent call last)" in text
                           and "Pipeline finished" not in text)
            if has_done and has_jobid:
                return ("SUCCESS", has_jobid.group(1))
            if has_crash:
                return ("FAILED", "exception before Pipeline finished")
        time.sleep(2)
    return ("TIMEOUT", f"no completion in {wait}s")

def _dump_failure(name, idx, workdir, log_file):
    if not log_file.exists():
        log.error("  [#%04d %s] pipeline.log missing.", idx, name); return
    leap = workdir / "build" / "leap.log"
    if leap.exists():
        text = leap.read_text()
        fatals = [ln for ln in text.splitlines() if "FATAL" in ln or "Error" in ln][-15:]
        if fatals:
            log.error("  [#%04d %s] leap.log FATAL/Error lines:", idx, name)
            for ln in fatals:
                log.error("    | %s", ln)
            return
    lines = [ln for ln in log_file.read_text().splitlines()
             if not ln.strip().startswith("[")
             and "Pre-flight" not in ln][-40:]
    log.error("  [#%04d %s] pipeline.log tail:", idx, name)
    for ln in lines:
        log.error("    | %s", ln)

def run_batch(args):
    batch_dir = Path(args.batch_dir).resolve()
    batch_dir.mkdir(parents=True, exist_ok=True)
    if not args.skip_preflight:
        log.info("Pre-flight: testing protein with tleap...")
        ok, msg = preflight_protein(args.protein_file, batch_dir)
        if not ok:
            log.error("="*70)
            log.error("PROTEIN PRE-FLIGHT FAILED:")
            for line in msg.splitlines():
                log.error("  %s", line)
            log.error("="*70)
            return 1
        log.info("  %s", msg)
    ligs = discover_ligands(args.ligands, batch_dir)
    if not ligs:
        log.error("No ligands found."); return 1
    log.info("Discovered %d ligand(s).", len(ligs))
    if args.dry_run:
        log.info("DRY RUN -- would process:")
        for idx, name, lig_file, charge in ligs:
            c = charge if charge is not None else args.ligand_charge
            log.info("  #%04d  %-40s  charge=%+d  %s", idx, name, c, lig_file.name)
        if args.decomp:
            mask = args.decomp_residues or "(all protein residues)"
            log.info("DECOMPOSITION enabled with mask: %s", mask)
        return 0
    user = os.environ.get("USER","")
    manifest, n_fail = [], 0
    for entry, (idx, name, lig_file, charge) in enumerate(ligs, 1):
        log.info("="*60)
        log.info("LIGAND %d/%d: #%04d %s", entry, len(ligs), idx, name)
        log.info("="*60)
        if args.max_concurrent > 0:
            _wait_for_capacity(user, args.queue, args.max_concurrent)
        proc, workdir, log_file = _submit_one(idx, name, lig_file, charge, args, batch_dir)
        status, detail = _parse_jobid(log_file, args.parse_jobid_timeout)
        if status == "SUCCESS":
            jobid = detail
            log.info("  -> LSF Job <%s>  pid=%d", jobid, proc.pid)
        else:
            jobid = status; n_fail += 1
            log.error("  -> %s (%s)  pid=%d", status, detail, proc.pid)
            _dump_failure(name, idx, workdir, log_file)
            if n_fail >= 5 and n_fail == entry:
                log.error("="*60)
                log.error("HALTING: first %d ligand(s) failed.", n_fail)
                log.error("="*60)
                _write_manifest(batch_dir, manifest +
                    [{"name":name,"ligfile":str(lig_file),"workdir":str(workdir),
                      "pid":proc.pid,"lsf_jobid":jobid,
                      "charge":charge if charge is not None else args.ligand_charge,
                      "submitted_at":datetime.now().isoformat()}])
                return 2
        manifest.append({"name":name,"ligfile":str(lig_file),
                         "workdir":str(workdir),"pid":proc.pid,"lsf_jobid":jobid,
                         "charge":charge if charge is not None else args.ligand_charge,
                         "submitted_at":datetime.now().isoformat()})
        _write_manifest(batch_dir, manifest)
    log.info("\nBatch complete: %d submitted (%d failed, %d ok).",
             len(manifest), n_fail, len(manifest)-n_fail)
    log.info("To aggregate: python -m amber_md.batch_aggregate %s", batch_dir)
    return 0 if n_fail == 0 else 1

def _write_manifest(batch_dir, manifest):
    tsv = batch_dir / "batch_manifest.tsv"
    with open(tsv, "w") as f:
        f.write("name\tlsf_jobid\tcharge\tworkdir\tligfile\tpid\tsubmitted_at\n")
        for m in manifest:
            f.write(f"{m['name']}\t{m['lsf_jobid']}\t{m['charge']}\t"
                    f"{m['workdir']}\t{m['ligfile']}\t{m['pid']}\t{m['submitted_at']}\n")
    (batch_dir / "batch_manifest.json").write_text(json.dumps(manifest, indent=2))

def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m amber_md.batch",
        description="Run the Amber MD pipeline over multiple ligands (v2.4.7).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--protein-file", type=Path, required=True)
    p.add_argument("--ligands", type=Path, required=True)
    p.add_argument("--batch-dir", type=Path, required=True)
    p.add_argument("--lig-resname", default="LIG")
    p.add_argument("--ligand-charge", type=int, default=0)
    p.add_argument("--charge-method", default="bcc", choices=["bcc","gas","resp"])
    p.add_argument("--salt", type=float, default=0.15)
    p.add_argument("--ion-method", default="rand", choices=["rand","grid"])
    p.add_argument("--no-gbsa", action="store_true")
    # v2.4.7: MM/GBSA per-residue decomposition (forwarded to each ligand's run_amber.py)
    p.add_argument("--decomp", action="store_true",
                   help="Enable MM/GBSA per-residue decomposition for every ligand "
                        "(produces FINAL_DECOMP_MMPBSA.dat per ligand).")
    p.add_argument("--decomp-residues", default="", metavar="MASK",
                   help="Amber residue mask for decomposition output, e.g. ':300-450'. "
                        "If empty with --decomp, MMPBSA defaults to all protein residues.")
    p.add_argument("--no-protonation", dest="protonation", action="store_false", default=True)
    p.add_argument("--prod-ns", type=float, default=10.0)  # final51: unified 10 ns
    p.add_argument("--equil-ns", type=float, default=1.0)
    p.add_argument("--project", default=None)
    p.add_argument("--queue", default="gpu")
    p.add_argument("--walltime", default="24:00")
    p.add_argument("--n-gpu", type=int, default=1)
    p.add_argument("--job-prefix", default="MD")
    p.add_argument("--max-concurrent", type=int, default=8)
    p.add_argument("--parse-jobid-timeout", type=int, default=300)
    p.add_argument("--skip-preflight", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--run-amber-py", type=Path, default=None)
    a = p.parse_args(argv)
    if a.run_amber_py is None:
        for c in [Path.cwd()/"run_amber.py",
                  Path(__file__).resolve().parent.parent/"run_amber.py"]:
            if c.exists(): a.run_amber_py = c; break
        else: p.error("Could not find run_amber.py.")
    log.info("="*60); log.info("AMBER MD BATCH SCREENING (v2.4.7)"); log.info("="*60)
    log.info("Protein:          %s", a.protein_file)
    log.info("Ligands:          %s", a.ligands)
    log.info("Batch dir:        %s", a.batch_dir)
    log.info("Concurrency cap:  %d", a.max_concurrent)
    log.info("Per-ligand:       %.1f ns equil + %.1f ns prod", a.equil_ns, a.prod_ns)
    log.info("LSF queue/proj:   %s / %s", a.queue, a.project or "(default)")
    log.info("RDKit available:  %s", "yes" if _have_rdkit() else "NO")
    if a.decomp:
        log.info("Decomposition:    ENABLED  mask=%s",
                 a.decomp_residues or "(all protein residues)")
    log.info("="*60)
    return run_batch(a)

if __name__ == "__main__":
    sys.exit(main() or 0)