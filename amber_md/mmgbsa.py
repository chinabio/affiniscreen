"""MMGBSAAnalyzer (v2.4.3) -- sane defaults for production-quality binding ΔG.

v2.4.3 changes:
  - Auto-generate an HTML report after a successful MM/GBSA run by calling
    mmpbsa_report.generate_report() on the FINAL_RESULTS_MMPBSA.dat file.
    Report generation is best-effort: any failure is logged as a warning and
    does NOT fail the run. The .dat file remains the contract.
  - New kwarg `make_report=True` (and CLI flag `--no-report`) to opt out,
    e.g. for minimal/headless environments without matplotlib.

v2.4.2 changes:
  - Default interval=2 (was 1): sample every other frame -- 2x speedup, no real
    accuracy loss for 10ns+ trajectories.
  - Default verbose=1 (was 2): per-frame summary instead of dumping all forces.
    This eliminates the reference.frc file (was hundreds of MB, hours of I/O).
  - Added netcdf=1: faster trajectory parsing.
  - Exposed --mmgbsa-interval, --mmgbsa-verbose, --mmgbsa-startframe via CLI.

Reference timings (10ns prod, ~75k atoms total):
  OLD defaults (interval=1, verbose=2):  ~2-3 hours, ~700MB intermediate files
  NEW defaults (interval=2, verbose=1):  ~10-15 minutes, ~50MB
  Report generation overhead:            ~1-2 seconds (negligible)
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import subprocess
from pathlib import Path
from .logger import get_logger
log = get_logger()


def write_mmgbsa_input(path: Path, *, igb: int = 8, saltcon: float = 0.15,
                       startframe: int = 1, endframe: int = 9999999,
                       interval: int = 2, verbose: int = 1,
                       decompose: bool = False) -> Path:
    """Generate a sane mmgbsa.in file. Returns the path written."""
    # IMPORTANT: MMPBSA.py's &general namelist parser is brittle about INLINE
    # comments -- a line like `interval=2,   ! sample every 2 frames` is parsed
    # as the value "2,! sample every 2 frames" and int() then raises
    # `invalid literal for int() with base 10`. So we emit a CLEAN namelist:
    # comments (if any) live on their OWN lines starting with `!`, never after
    # a value on the same line.
    body = f"""MM/GBSA binding free energy (v2.4.2 defaults)
! interval = sample every N frames; verbose 1 = per-frame, 2 = + per-atom (slow)
! igb 8 = GBn2 (Onufriev/Bashford/Case 2); saltcon in M (physiological ~0.15)
&general
  startframe={startframe},
  endframe={endframe},
  interval={interval},
  verbose={verbose},
  keep_files=1,
  netcdf=1,
/
&gb
  igb={igb},
  saltcon={saltcon},
/
"""
    if decompose:
        body += """&decomp
  idecomp=2,
  dec_verbose=1,
  print_res="within 5",
/
"""
    path.write_text(body)
    return path


class MMGBSAAnalyzer:
    """Run MMPBSA.py on a trajectory using the existing per-component topologies."""

    def __init__(self, topo_dir, traj_file, workdir, *,
                 igb=8, saltcon=0.15, interval=2, verbose=1,
                 startframe=1, decompose=False, make_report=True):
        self.topo_dir = Path(topo_dir).resolve()
        self.traj_file = Path(traj_file).resolve()
        self.workdir = Path(workdir).resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.igb = igb
        self.saltcon = saltcon
        self.interval = interval
        self.verbose = verbose
        self.startframe = startframe
        self.decompose = decompose
        self.make_report = make_report

    def run(self):
        log.info("MMGBSAAnalyzer: using existing topologies in %s", self.topo_dir)
        sp = self.topo_dir / "solvated.prmtop"
        cp = self.topo_dir / "complex.prmtop"
        rp = self.topo_dir / "receptor.prmtop"
        lp = self.topo_dir / "ligand.prmtop"
        for f in (sp, cp, rp, lp, self.traj_file):
            if not f.exists():
                raise FileNotFoundError(f"Required input missing: {f}")

        infile = self.workdir / "mmgbsa.in"
        write_mmgbsa_input(infile, igb=self.igb, saltcon=self.saltcon,
                           startframe=self.startframe, interval=self.interval,
                           verbose=self.verbose, decompose=self.decompose)

        outfile = self.workdir / "FINAL_RESULTS_MMPBSA.dat"
        cmd = ["MMPBSA.py", "-O",
               "-i", str(infile),
               "-o", str(outfile),
               "-sp", str(sp), "-cp", str(cp),
               "-rp", str(rp), "-lp", str(lp),
               "-y", str(self.traj_file)]
        if self.decompose:
            cmd += ["-do", str(self.workdir / "decomp_per_res.dat")]

        log.info("RUN: %s  (cwd=%s)", " ".join(cmd), self.workdir)
        cp_res = subprocess.run(cmd, cwd=str(self.workdir),
                                capture_output=True, text=True)
        if cp_res.returncode != 0:
            log.error("MMPBSA.py failed (rc=%d). stderr:\n%s",
                      cp_res.returncode, cp_res.stderr[-2000:])
            # Self-diagnose: MMPBSA.py leaves per-calculation temp files whose
            # tails carry the REAL error (e.g. the sander/mmpbsa_py_energy
            # mdout). Surface them so we don't have to guess from a generic
            # CalcError. keep_files=0 still leaves them on a failed run.
            try:
                import glob as _glob, shutil as _shutil
                _keep = self.workdir / "mmpbsa_failed_tmp"
                try: _keep.mkdir(exist_ok=True)
                except Exception: pass
                diag = sorted(
                    _glob.glob(str(self.workdir / "_MMPBSA_*.out*")) +
                    _glob.glob(str(self.workdir / "_MMPBSA_*.mdout*")) +
                    _glob.glob(str(self.workdir / "reference.frc")))
                # Preserve a copy before MPI/cleanup can delete them.
                for fp in diag:
                    try: _shutil.copy2(fp, _keep / Path(fp).name)
                    except Exception: pass
                if diag:
                    log.error("Retained MMPBSA temp files in: %s", _keep)
                for fp in diag[:6]:
                    try:
                        tail = "".join(open(fp).read().splitlines(keepends=True)[-25:])
                        log.error("---- %s (tail) ----\n%s", Path(fp).name, tail)
                    except Exception:  # noqa: BLE001
                        pass
                if not diag:
                    log.error("No _MMPBSA_*.out temp files found in %s to "
                              "diagnose; check radii/igb consistency.",
                              self.workdir)
            except Exception as _e:  # noqa: BLE001
                log.error("Could not collect MMPBSA temp-file diagnostics: %s", _e)
            raise RuntimeError(f"MMPBSA.py failed with rc={cp_res.returncode}")
        if not outfile.exists():
            raise RuntimeError(f"MMPBSA.py finished but no output at {outfile}")

        # Quick parse of the bottom line for the log
        text = outfile.read_text()
        import re
        m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", text)
        if m:
            log.info("ΔG_bind = %s ± %s kcal/mol (SEM %s)",
                     m.group(1), m.group(2), m.group(3))

        # Best-effort HTML report. Never fails the MM/GBSA run.
        if self.make_report:
            self._write_report(outfile)

        return outfile

    def _write_report(self, dat_path):
        """Generate a styled HTML report next to the .dat file.

        Best-effort: any failure (missing matplotlib, parse error, write error)
        is logged as a warning and swallowed. The .dat file is the contract;
        the HTML is a convenience.
        """
        try:
            # Lazy import so MMGBSAAnalyzer stays usable on minimal/headless
            # environments that don't have matplotlib installed.
            from .mmpbsa_report import generate_report
            html_path = generate_report(dat_path)
            log.info("Report written: %s", html_path)
        except Exception as e:
            log.warning("Report generation skipped (%s: %s)",
                        type(e).__name__, e)
