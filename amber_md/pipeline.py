"""End-to-end orchestrator (v2.5.1 - in-package analysis_kit / mmpbsa_report)."""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .config import WorkflowConfig
from .logger import get_logger
from .utils import ensure_dir, CommandError
from .prep import PDBCleaner, LigandParametrizer
from .protonation import apply_protonation
from .builder import SystemBuilder
from .md_inputs import MDInputWriter
from .submit import LSFSubmitter
from .monitor import JobMonitor
from .analysis import TrajectoryAnalyzer, MMGBSAAnalyzer
from .topology import TopologySplitter
from .visualize import Visualizer


class AmberPipeline:
    def __init__(self, cfg):
        self.cfg = cfg
        self.work = ensure_dir(Path(cfg.work_dir).resolve())
        self.log  = get_logger(log_file=self.work/"pipeline.log")
        cfg.save(self.work/"workflow_config.json")

    def stage_prep(self):
        self.log.info("=== Stage 1/5: PREP ===")
        c = PDBCleaner(self.work/"prep")

        if self.cfg.protein_input and self.cfg.ligand_input:
            self.log.info("Mode: DIRECT (separate protein + ligand files)")
            self.log.info("  protein: %s", self.cfg.protein_input)
            self.log.info("  ligand : %s", self.cfg.ligand_input)
            prot = c.clean_protein_only(self.cfg.protein_input)
            if getattr(self.cfg, "auto_protonation", True):
                prot_fixed = prot.parent / (prot.stem + "_proton.pdb")
                apply_protonation(prot, prot_fixed,
                    manual_overrides=getattr(self.cfg, "protonation_overrides", None))
                prot = prot_fixed
            lig_pdb = Path(self.cfg.ligand_input)
        else:
            self.log.info("Mode: COMBINED (single complex PDB, split by resname)")
            if not self.cfg.complex_pdb or not Path(self.cfg.complex_pdb).exists():
                raise RuntimeError(
                    "COMBINED mode requires --pdb pointing at a real file, "
                    "OR provide --protein-file and --ligand-file for DIRECT mode.")
            clean = c.clean(self.cfg.complex_pdb)
            prot, lig_pdb = c.split_complex(
                clean, self.cfg.ligand_resname, self.work/"prep")
            if getattr(self.cfg, "auto_protonation", True):
                prot_fixed = prot.parent / (prot.stem + "_proton.pdb")
                apply_protonation(prot, prot_fixed,
                    manual_overrides=getattr(self.cfg, "protonation_overrides", None))
                prot = prot_fixed
            if self.cfg.ligand_input:
                lig_pdb = Path(self.cfg.ligand_input)

        lp = LigandParametrizer(self.work, self.cfg.system)
        mol2, frc = lp.parametrize(lig_pdb, self.cfg.ligand_resname)
        return prot, mol2, frc

    def stage_build(self, prot, lig_mol2, lig_frcmod):
        self.log.info("=== Stage 2/5: BUILD ===")
        builder = SystemBuilder(self.work, self.cfg.system)
        prmtop, inpcrd, pdb = builder.build(
            prot, lig_mol2, lig_frcmod, self.cfg.ligand_resname)
        try:
            sp = TopologySplitter(self.work, self.cfg.ligand_resname)
            sp.sanity_check(sp.split(prmtop))
        except Exception as e:
            self.log.warning("Topology split skipped: %s", e)
        try:
            Visualizer(self.work).render_all(pdb, ligand_resname=self.cfg.ligand_resname)
        except Exception as e:
            self.log.warning("Pre-MD viz skipped: %s", e)
        return prmtop, inpcrd

    def stage_md(self, prmtop, inpcrd):
        self.log.info("=== Stage 3/5: MD ===")
        writer = MDInputWriter(self.work, self.cfg.md)
        mdin = writer.write_all()

        if self.cfg.mmgbsa.enabled:
            try:
                ana = MMGBSAAnalyzer(
                    self.work, self.cfg.mmgbsa,
                    ligand_resname=self.cfg.ligand_resname)
                ana.prepare_inputs()
            except Exception as e:
                self.log.warning(
                    "MM/GBSA input prep failed (in-job MM/GBSA will be skipped): %s", e)

        subm = LSFSubmitter(self.work/"jobs", self.cfg.hpc)
        script = subm.build_gpu_md_script(
            mdin, prmtop, inpcrd,
            mmgbsa_cfg=self.cfg.mmgbsa,
            ligand_resname=self.cfg.ligand_resname,
            workdir=self.work,
        )

        if not self.cfg.submit:
            self.log.info("submit=False -> manual launch.")
            return None
        jid = subm.submit(script)
        if self.cfg.monitor:
            JobMonitor().wait(jid)
        return jid

    def stage_analyze(self, prmtop, job_id=None):
        self.log.info("=== Stage 4/5: ANALYSIS ===")
        traj = self.work/"jobs"/"prod.nc"
        if not traj.exists():
            if not self.cfg.submit:
                self.log.info("Analysis skipped: MD was not submitted (submit=False).")
                self._print_resume_hint(job_id=None)
            elif not self.cfg.monitor:
                self.log.info(
                    "Analysis skipped: GPU job is still running (--no-monitor was set).")
                self._print_resume_hint(job_id=job_id)
            else:
                self.log.warning("Trajectory %s not found (job may have failed).", traj)
                self.log.warning("Check job logs in %s for errors.", self.work/"jobs")
                self._print_resume_hint(job_id=job_id)
            return

        ta = TrajectoryAnalyzer(self.work, prmtop, traj)
        rmsd, rmsf = ta.rmsd_rmsf()
        ta.plot(rmsd, rmsf)
        try:
            Visualizer(self.work).render_all(
                prmtop, traj, ligand_resname=self.cfg.ligand_resname)
        except Exception as e:
            self.log.warning("Post-MD viz skipped: %s", e)

        if self.cfg.mmgbsa.enabled:
            final = self.work/"mmgbsa"/"FINAL_RESULTS_MMPBSA.dat"
            if final.exists() and final.stat().st_size > 0:
                self.log.info("MM/GBSA already complete (from LSF job): %s", final)
            else:
                self.log.info(
                    "MM/GBSA result not found from LSF job - running locally as fallback.")
                ana = MMGBSAAnalyzer(
                    self.work, self.cfg.mmgbsa,
                    ligand_resname=self.cfg.ligand_resname)
                try:
                    ana.run_serial(solvated_prm=prmtop, traj=traj)
                except CommandError as e:
                    self.log.error("MM/GBSA failed: %s", e)

            self._ensure_mmgbsa_report(final)

        self._ensure_analysis_kit()

    def _ensure_mmgbsa_report(self, dat_path):
        """v2.5.1: prefers package-relative import, falls back to file discovery
        first inside the package, then at the legacy sibling location.
        Never raises - failures are logged but never fail the pipeline."""
        if not dat_path.exists() or dat_path.stat().st_size == 0:
            return
        html_path = dat_path.with_name("FINAL_RESULTS.report.html")
        if html_path.exists() and html_path.stat().st_size > 0:
            self.log.info("HTML report already present: %s", html_path)
            return
        try:
            generate_report = self._load_generate_report()
            if generate_report is None:
                self.log.info("mmpbsa_report not importable - skipping HTML report.")
                return
            out = generate_report(dat_path, html_path)
            self.log.info("Wrote HTML report: %s", out)
        except Exception as e:
            self.log.warning("HTML report generation failed (non-fatal): %s", e)

    def _load_generate_report(self):
        """Resolve mmpbsa_report.generate_report.
        Order: package import -> in-package file -> legacy sibling file."""
        try:
            from .mmpbsa_report import generate_report  # type: ignore
            return generate_report
        except Exception:
            pass

        import importlib.util, sys as _sys
        pkg_dir = Path(__file__).resolve().parent
        for cand in (pkg_dir / "mmpbsa_report.py",
                     pkg_dir.parent / "mmpbsa_report.py"):
            if not cand.exists():
                continue
            try:
                spec = importlib.util.spec_from_file_location("mmpbsa_report", cand)
                mod = importlib.util.module_from_spec(spec)
                _sys.modules["mmpbsa_report"] = mod
                spec.loader.exec_module(mod)
                return getattr(mod, "generate_report", None)
            except Exception:
                continue
        return None

    def _ensure_analysis_kit(self):
        """v2.5.1: run analysis_kit/run_analysis.sh if outputs are missing.
        In-package first, then legacy fallback. Non-fatal."""
        report = self.work / "analysis" / "COMBINED_REPORT.html"
        if report.exists() and report.stat().st_size > 0:
            self.log.info("Analysis kit output already present: %s", report)
            return

        ak = self._find_analysis_kit()
        if ak is None:
            self.log.info("analysis_kit not found - skipping kit run.")
            return
        runner = ak / "run_analysis.sh"
        if not runner.exists():
            self.log.info("analysis_kit/run_analysis.sh not found at %s.", runner)
            return

        prod = self.work / "jobs" / "prod.nc"
        cprm = self.work / "build" / "complex.prmtop"
        if not (prod.exists() and cprm.exists()):
            self.log.info("Skipping analysis_kit: prod.nc or complex.prmtop missing.")
            return

        import subprocess
        log_file = self.work / "analysis_kit.log"
        self.log.info("Running analysis_kit on %s (log: %s)", self.work, log_file)
        try:
            with open(log_file, "w") as f:
                cp = subprocess.run(
                    ["bash", str(runner), str(self.work)],
                    cwd=str(ak), stdout=f, stderr=subprocess.STDOUT, timeout=900)
            if cp.returncode == 0:
                self.log.info("analysis_kit done -> %s", report)
            else:
                self.log.warning("analysis_kit exited %d (see %s)",
                                 cp.returncode, log_file)
        except subprocess.TimeoutExpired:
            self.log.warning("analysis_kit timed out (900 s) - see %s", log_file)
        except Exception as e:
            self.log.warning("analysis_kit invocation failed (non-fatal): %s", e)

    def _find_analysis_kit(self):
        """In-package first, then legacy sibling fallback."""
        pkg_dir = Path(__file__).resolve().parent
        for cand in (pkg_dir / "analysis_kit", pkg_dir.parent / "analysis_kit"):
            if cand.is_dir():
                return cand
        return None

    def _print_resume_hint(self, job_id=None):
        workdir = self.work
        msg = ["", "  " + "="*64,
               "  TO RUN ANALYSIS WHEN THE JOB COMPLETES:",
               "  " + "="*64]
        if job_id:
            msg.append(f"  1. Wait for job {job_id} to finish:")
            msg.append(f"       bjobs {job_id}        # status (RUN -> DONE)")
            msg.append(f'       bwait -w "done({job_id})"   # block until done')
        else:
            msg.append("  1. Wait for the GPU job to finish:")
            msg.append("       bjobs            # check status")
        msg.append("")
        msg.append("  2. Verify the trajectory exists:")
        msg.append(f"       ls -lh {workdir}/jobs/prod.nc")
        msg.append("")
        msg.append("  3. Re-run just the analysis stage:")
        msg.append(f"       python -m amber_md.resume {workdir}")
        msg.append("")
        msg.append("     (this skips PREP/BUILD/MD, runs Stage 4 + 5 only)")
        msg.append("  " + "="*64)
        for line in msg:
            self.log.info(line)

    def stage_fep(self):
        if not getattr(self.cfg, "fep", None) or not self.cfg.fep.enabled:
            return
        self.log.info("=== Stage 5/5: ALCHEMICAL FEP ===")
        from .fep import FEPSetup, FEPAnalyzer, relative_binding_dG
        fep = FEPSetup(self.work, self.cfg.fep.params, self.cfg.md, self.cfg.hpc)
        legs = {}
        if self.cfg.fep.complex_prmtop and self.cfg.fep.complex_inpcrd:
            legs["complex"] = fep.setup_leg(
                "complex",
                Path(self.cfg.fep.complex_prmtop), Path(self.cfg.fep.complex_inpcrd))
        if self.cfg.fep.solvent_prmtop and self.cfg.fep.solvent_inpcrd:
            legs["solvent"] = fep.setup_leg(
                "solvent",
                Path(self.cfg.fep.solvent_prmtop), Path(self.cfg.fep.solvent_inpcrd))
        jids = {}
        if self.cfg.submit:
            for n, d in legs.items():
                jids[n] = fep.submit_leg(d, n)
            if self.cfg.monitor:
                mon = JobMonitor()
                for jid in jids.values():
                    mon.wait(jid)
        results = {n: FEPAnalyzer(d, self.cfg.fep.params.lambdas).run()
                   for n, d in legs.items()}
        if "complex" in results and "solvent" in results:
            ddG = relative_binding_dG(results["complex"], results["solvent"])
            self.log.info("ddG_bind = %s kcal/mol",
                          f"{ddG:.3f}" if ddG is not None else "N/A")
            (self.work/"fep"/"RESULTS.txt").write_text(
                f"dG_complex = {results['complex']['dG_kcal_mol']}\n"
                f"dG_solvent = {results['solvent']['dG_kcal_mol']}\n"
                f"ddG_bind   = {ddG}\n")
        return results

    def run(self):
        try:
            prot, mol2, frc = self.stage_prep()
            prm, crd        = self.stage_build(prot, mol2, frc)
            jid             = self.stage_md(prm, crd)
            self.stage_analyze(prm, job_id=jid)
            self.stage_fep()
            self.log.info("Pipeline finished.")
        except CommandError as e:
            self.log.exception("Pipeline failed: %s", e)
            raise
