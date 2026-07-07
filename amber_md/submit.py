"""LSF submission (v2.5.3).

v2.5.3: _header() emits `set -euo pipefail` BEFORE module load /
source, and resolves venv_activate to an absolute path with an
existence guard. Fixes silent activate-script failure that left
MM/GBSA without its environment.

v2.5.2: export PYTHONPATH=<workflow_root> in the LSF tail block before the
in-job `python -c "from amber_md.mmpbsa_report import ..."` call. Without
this, the compute node's Python could not find `amber_md` because cwd on
the node is wherever LSF drops the job, not the workflow root, and the
activate script does not put the workflow root on PYTHONPATH itself.
Symptom: tail block printed
    ModuleNotFoundError: No module named 'amber_md'
and the per-ligand HTML report was silently never generated.

v2.5.1: _analysis_tail_block resolves analysis_kit INSIDE the amber_md
package (Path(__file__).parent / "analysis_kit") where it actually lives,
with a fallback to the legacy parent.parent location.

v2.5.0:
  * _analysis_tail_block() appended to the GPU MD script after MM/GBSA.
  * build_gpu_md_script gained 'run_analysis' kwarg (default True).

v2.4.9:
  * Report generation calls the package API directly.
  * `set +e` / `set -e` wrapper around the report call.

v2.4.8:
  * In-job MM/GBSA tail-block also generates an HTML report.
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .config import HPCConfig
from .utils import run, ensure_dir
from .logger import get_logger
log = get_logger()

GPU_MEM_MB_PER_SLOT = 4096

# Package name for in-LSF `python -c "from <pkg>.mmpbsa_report import ..."`
_PKG_NAME = __package__ or "amber_md"

# v2.5.2: workflow root = directory containing the amber_md/ package.
# Exported as PYTHONPATH in the LSF tail block so `import amber_md` works
# on compute nodes where cwd is not the workflow root.
_WORKFLOW_ROOT = str(Path(__file__).resolve().parent.parent)


def _default_analysis_kit_dir():
    """v2.5.1: in-package layout first, then legacy sibling fallback."""
    pkg_dir = Path(__file__).resolve().parent
    in_pkg = pkg_dir / "analysis_kit"
    if in_pkg.is_dir():
        return in_pkg
    legacy = pkg_dir.parent / "analysis_kit"
    if legacy.is_dir():
        return legacy
    return in_pkg


class LSFSubmitter:
    def __init__(self, work_dir, cfg):
        self.work_dir = ensure_dir(work_dir)
        self.cfg = cfg

    def _header(self, queue, gpu, name, log_prefix, n_slots=None):
        h  = "#!/bin/bash\n"
        h += f"#BSUB -q {queue}\n#BSUB -P {self.cfg.project}\n"
        h += f"#BSUB -J {name}\n#BSUB -W {self.cfg.walltime}\n"
        h += f"#BSUB -o {log_prefix}.%J.out\n#BSUB -e {log_prefix}.%J.err\n"
        if gpu:
            ngpu = int(getattr(self.cfg, "n_gpu", 1) or 1)
            mode = str(getattr(self.cfg, "gpu_request_mode", "rusage")).lower()
            mem = GPU_MEM_MB_PER_SLOT
            if n_slots is not None:
                # Explicit caller override (e.g. FEP). Honor it verbatim with
                # -gpu num=ngpu; caller is responsible for slot semantics.
                h += f"#BSUB -n {int(n_slots)}\n"
                h += '#BSUB -R "span[hosts=1]"\n'
                h += f'#BSUB -R "rusage[mem={mem}]"\n'
                h += f'#BSUB -gpu "num={ngpu}"\n'
            elif mode == "rusage":
                # final54: slots==GPUs queues (the cluster `gpu`). Request CPU cores
                # via -n, and the GPU(s) as a RESOURCE, so cores never inflate
                # the GPU count. Enables parallel MM/GBSA scoring on 1 GPU.
                cores = int(getattr(self.cfg, "gpu_cpu_cores", None) or 1)
                h += f"#BSUB -n {cores}\n"
                h += '#BSUB -R "span[hosts=1]"\n'
                h += f'#BSUB -R "rusage[mem={mem},ngpus_physical={ngpu}]"\n'
                # final55: this LSF rejects GPU in BOTH rusage and -gpu; request
                # the GPU ONLY via rusage[ngpus_physical] -- no -gpu line here.'
            else:  # mode == "slots" (legacy: slots==cores)
                slots = int(getattr(self.cfg, "n_gpu_slots", None) or ngpu or 1)
                h += f"#BSUB -n {slots}\n"
                h += '#BSUB -R "span[hosts=1]"\n'
                h += f'#BSUB -R "rusage[mem={mem}]"\n'
                h += f'#BSUB -gpu "num={ngpu}"\n'
        else:
            h += f"#BSUB -n {self.cfg.n_cpu}\n#BSUB -R span[ptile={self.cfg.n_cpu}]\n"
        # v2.5.3 FIX: enable strict mode BEFORE module load / source so a
        # missing activate script (or any failing step) aborts the job
        # immediately instead of silently continuing into a doomed
        # MM/GBSA stage.
        h += "\nset -euo pipefail\n"
        h += "module purge\n"
        for m in self.cfg.modules:
            h += f"module load {m}\n"
        if self.cfg.venv_activate:
            # v2.5.3 FIX: resolve the activate script to an ABSOLUTE path.
            # LSF runs the job from ~/.lsbatch/..., not the workflow root,
            # so a relative "./activate_amber_md.sh" is never found.
            # Resolving here covers EVERY launch path (GUI/batch/run_amber).
            _act = Path(self.cfg.venv_activate)
            if not _act.is_absolute():
                _cand = Path(_WORKFLOW_ROOT) / _act.name
                _act = _cand if _cand.exists() else _act
            h += f'if [ ! -f "{_act}" ]; then\n'
            h += f'    echo "ERROR: activate script not found: {_act}" >&2\n'
            h += '    exit 1\n'
            h += 'fi\n'
            h += f"source {_act}\n"
        # v2.5.2: ensure the workflow root is on PYTHONPATH for the whole
        # job, not just the tail block. Belt-and-braces: the tail block
        # also re-exports it locally.
        h += f'export PYTHONPATH="{_WORKFLOW_ROOT}${{PYTHONPATH:+:$PYTHONPATH}}"\n\n'
        return h

    def build_gpu_md_script(self, mdin, prmtop, inpcrd, script_path=None,
                            mmgbsa_cfg=None, ligand_resname="LIG", workdir=None,
                            run_analysis=True):
        sp = script_path or self.work_dir/"run_md.lsf"
        b = self._header(self.cfg.queue_gpu, True, self.cfg.job_name, "md")
        b += f"PRM={prmtop}\nCRD={inpcrd}\n\n"
        b += f"pmemd.cuda -O -i {mdin['min1']} -p $PRM -c $CRD -ref $CRD -o min1.out -r min1.rst -x min1.nc\n"
        b += f"pmemd.cuda -O -i {mdin['min2']} -p $PRM -c min1.rst -o min2.out -r min2.rst -x min2.nc\n"
        b += f"pmemd.cuda -O -i {mdin['heat']} -p $PRM -c min2.rst -ref min2.rst -o heat.out -r heat.rst -x heat.nc\n"
        b += f"pmemd.cuda -O -i {mdin['equil']} -p $PRM -c heat.rst -ref heat.rst -o equil.out -r equil.rst -x equil.nc\n"
        b += f"pmemd.cuda -O -i {mdin['prod']} -p $PRM -c equil.rst -o prod.out -r prod.rst -x prod.nc\n"

        if mmgbsa_cfg is not None and getattr(mmgbsa_cfg, "enabled", False):
            b += self._mmgbsa_tail_block(mmgbsa_cfg, ligand_resname, workdir)

        if run_analysis and workdir is not None:
            b += self._analysis_tail_block(workdir)

        sp.write_text(b)
        sp.chmod(0o755)
        return sp

    def _mmgbsa_tail_block(self, mmgbsa_cfg, ligand_resname, workdir):
        if workdir is None:
            raise ValueError(
                "build_gpu_md_script: 'workdir' is required when MM/GBSA is enabled.")
        wd = str(Path(workdir).resolve())
        cap = int(getattr(self.cfg, "mmgbsa_n_cpu", None) or 0)

        decomp_arg = ""
        if getattr(mmgbsa_cfg, "decomposition", False):
            decomp_arg = f" -do {wd}/mmgbsa/FINAL_DECOMP_MMPBSA.dat"

        block = [
            "",
            "# ------------------------------------------------------------------",
            "# v2.4.9: MM/GBSA on the GPU node's CPU cores, right after MD.",
            "# v2.5.2: PYTHONPATH is set in the header so `import amber_md`",
            "# works on the compute node (cwd != workflow root). The local",
            "# export below is belt-and-braces in case the header changes.",
            "# ------------------------------------------------------------------",
            f"cd {wd}/jobs",
            "if [ ! -s prod.nc ]; then",
            '    echo "[MM/GBSA] prod.nc missing or empty - skipping." >&2',
            "    exit 0",
            "fi",
            "",
            'NP="${LSB_DJOB_NUMPROC:-1}"',
            f"MMGBSA_CAP={cap}",
            'if [ "$MMGBSA_CAP" -gt 0 ] && [ "$MMGBSA_CAP" -lt "$NP" ]; then',
            '    NP="$MMGBSA_CAP"',
            "fi",
            'if [ "$NP" -gt 1 ]; then',
            '    MMPBSA_CMD="mpirun -np ${NP} MMPBSA.py.MPI"',
            "else",
            '    MMPBSA_CMD="MMPBSA.py"',
            "fi",
            'echo "[MM/GBSA] using ${MMPBSA_CMD}  (LSB_DJOB_NUMPROC=${LSB_DJOB_NUMPROC:-unset}, cap=${MMGBSA_CAP})"',
            "",
            f"mkdir -p {wd}/mmgbsa",
            f"cd {wd}/mmgbsa",
            "",
            "${MMPBSA_CMD} -O \\",
            f"  -i  {wd}/mmgbsa/mmgbsa.in \\",
            f"  -o  {wd}/mmgbsa/FINAL_RESULTS_MMPBSA.dat \\",
            f"  -sp {wd}/build/complex.prmtop \\",
            f"  -cp {wd}/topo/complex.prmtop \\",
            f"  -rp {wd}/topo/receptor.prmtop \\",
            f"  -lp {wd}/topo/ligand.prmtop \\",
            f"  -y  {wd}/jobs/prod.nc{decomp_arg} \\",
            f"  > {wd}/mmgbsa/mmgbsa.log 2>&1",
            "",
            f'echo "[MM/GBSA] done $(date)" >> {wd}/mmgbsa/mmgbsa.log',
            "",
            f"if [ -s {wd}/mmgbsa/FINAL_RESULTS_MMPBSA.dat ]; then",
            '    echo "[Report] generating HTML report..."',
            "    set +e",
            # v2.5.2: re-export PYTHONPATH locally (header already does it,
            # but this makes the tail block self-contained for resumes).
            f'    export PYTHONPATH="{_WORKFLOW_ROOT}${{PYTHONPATH:+:$PYTHONPATH}}"',
            f'    python -c "from {_PKG_NAME}.mmpbsa_report import generate_report; "\\',
            f'"generate_report(\'{wd}/mmgbsa/FINAL_RESULTS_MMPBSA.dat\')" \\',
            f"        >> {wd}/mmgbsa/mmgbsa.log 2>&1",
            "    rc=$?",
            "    set -e",
            "    if [ $rc -eq 0 ]; then",
            f'        echo "[Report] wrote {wd}/mmgbsa/FINAL_RESULTS.report.html"',
            "    else",
            f'        echo "[Report] generation failed (rc=$rc, non-fatal) - see {wd}/mmgbsa/mmgbsa.log"',
            "    fi",
            "fi",
            "",
        ]
        return "\n".join(block)

    def _analysis_tail_block(self, workdir, analysis_kit_dir=None):
        """v2.5.1: in-package analysis_kit by default."""
        wd = str(Path(workdir).resolve())
        if analysis_kit_dir is None:
            analysis_kit_dir = _default_analysis_kit_dir()
        ak = str(Path(analysis_kit_dir).resolve())

        block = [
            "",
            "# ------------------------------------------------------------------",
            "# v2.5.0: post-MD analysis kit (cpptraj RMSD/RMSF/hbond/contacts +",
            "# PyMOL/VMD load scripts + COMBINED_REPORT.html).",
            "# Non-fatal: a failure here will NOT fail the LSF job.",
            "# ------------------------------------------------------------------",
            f"if [ -s {wd}/jobs/prod.nc ] && [ -s {wd}/build/complex.prmtop ]; then",
            f"    if [ -x {ak}/run_analysis.sh ]; then",
            '        echo "[Analysis] launching analysis kit..."',
            "        set +e",
            f"        bash {ak}/run_analysis.sh {wd} > {wd}/analysis_kit.log 2>&1",
            "        rc=$?",
            "        set -e",
            "        if [ $rc -eq 0 ]; then",
            f'            echo "[Analysis] done -> {wd}/analysis/COMBINED_REPORT.html"',
            "        else",
            f'            echo "[Analysis] failed (rc=$rc, non-fatal) - see {wd}/analysis_kit.log"',
            "        fi",
            "    else",
            f'        echo "[Analysis] skipped: {ak}/run_analysis.sh not found / not executable" >&2',
            "    fi",
            "else",
            '    echo "[Analysis] skipped: prod.nc or complex.prmtop missing" >&2',
            "fi",
            "",
        ]
        return "\n".join(block)

    def build_cpu_mmgbsa_script(self, command, script_path=None):
        sp = script_path or self.work_dir/"run_mmgbsa.lsf"
        b = self._header(self.cfg.queue_cpu, False,
                         self.cfg.job_name+"_gbsa", "mmgbsa")
        b += command + "\n"
        sp.write_text(b)
        sp.chmod(0o755)
        return sp

    def submit(self, script):
        return _bsub_submit(
            Path(script), self.work_dir,
            project=self.cfg.project,
            queue=(self.cfg.queue_gpu if "md" in Path(script).name.lower()
                   else self.cfg.queue_cpu),
            walltime=self.cfg.walltime)


def _bsub_submit(script, cwd, *, project, queue, walltime, extra_args=None):
    import subprocess
    cmd = ["bsub", "-P", project, "-q", queue, "-W", walltime]
    if extra_args:
        cmd.extend(extra_args)
    log.info("RUN: %s (stdin=%s)  (cwd=%s)", " ".join(cmd), script, cwd)
    script_text = Path(script).read_text()
    try:
        cp = subprocess.run(cmd, cwd=str(cwd), input=script_text,
                            text=True, capture_output=True, check=False)
    except FileNotFoundError as e:
        raise RuntimeError(f"bsub not found: {e}") from e
    if cp.returncode != 0:
        raise RuntimeError(
            f"bsub failed (exit {cp.returncode}):\n"
            f"STDOUT: {cp.stdout}\nSTDERR: {cp.stderr}\n"
            f"Script: {script}")
    log.info("bsub stdout: %s", (cp.stdout or "").strip())
    import re
    m = re.search(r"Job <(\d+)>", cp.stdout or "")
    if not m:
        raise RuntimeError(
            f"Could not parse bsub output (expected 'Job <NNN>'):\n{cp.stdout}")
    return m.group(1)
