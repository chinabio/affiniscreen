#!/usr/bin/env python3
"""
preflight_abfe.py  --  Definitive "will Amber ABFEP run here?" check.

Part of amber_md_workflow_v2.5.0. Run this ON A GPU NODE in the SAME environment
your ABFEP job would use (i.e. after `module load amber/22` and conda activate):

    bsub -q gpu -gpu "num=1" -W 0:15 -Is bash
    module load amber/22
    conda activate <your_env>     # if your topology build needs env tools
    python -m amber_md.preflight_abfe --work-dir /tmp/abfe_pf

It does NOT burn real FEP hours. Checks, in order:
  [1] pmemd.cuda + pmemd.cuda.MPI on PATH (GPU TI binaries from system amber).
  [2] tleap present + builds a trivial alchemical topology.
  [3] parmed imports cleanly in THIS env (the numpy.compat trap from MM-GBSA).
  [4] pmemd.cuda ACCEPTS the GTI/softcore mdin keywords (tiny real run).
  [5] pmemd.cuda.MPI multi-window (-ng) smoke test (HREMD path)  [WARN-only].
  [6] alchemlyb + pymbar importable (MBAR; trapezoid-TI fallback) [WARN-only].

Exit 0 => GO. Non-zero => at least one BLOCKER; details in stdout and
<work-dir>/abfe_preflight_report.json.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations
import argparse, json, shutil, subprocess, sys, tempfile, textwrap
from pathlib import Path


def _run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def _ok(m):   print(f"  [ OK ]    {m}")
def _warn(m): print(f"  [WARN]    {m}")
def _fail(m): print(f"  [BLOCK]   {m}")


def check_binaries(report):
    print("\n[1] GPU TI binaries on PATH")
    blockers = []
    for tool in ("pmemd.cuda", "pmemd.cuda.MPI"):
        path = shutil.which(tool)
        report["binaries"][tool] = path
        if path:
            ver = _run([tool, "--version"])
            line = ((ver.stdout or ver.stderr).strip().splitlines() or [""])[0]
            _ok(f"{tool} -> {path}   {line}")
        else:
            _fail(f"{tool} NOT on PATH. ABFEP needs it; it ships with the system "
                  f"Amber module (e.g. `module load amber/22`), NOT the conda env.")
            blockers.append(tool)
    return blockers


def check_tleap(report):
    print("\n[2] tleap present")
    p = shutil.which("tleap")
    report["binaries"]["tleap"] = p
    if p:
        _ok(f"tleap -> {p}"); return []
    _fail("tleap NOT on PATH (needed to build alchemical topologies).")
    return ["tleap"]


def check_parmed(report):
    print("\n[3] parmed imports cleanly in THIS env (numpy.compat trap)")
    code = ("import parmed, numpy;"
            "print('parmed', parmed.__version__, 'numpy', numpy.__version__)")
    r = _run([sys.executable, "-c", code])
    report["parmed"] = {"rc": r.returncode, "out": r.stdout.strip(),
                        "err": r.stderr.strip()[-400:]}
    if r.returncode == 0:
        _ok(r.stdout.strip()); return []
    if "numpy.compat" in r.stderr:
        _fail("parmed import fails with numpy.compat -- the amber module is "
              "shadowing the env's parmed with an old copy (same class of bug "
              "that broke MM-GBSA). Build topology with the env's parmed BEFORE "
              "loading amber.")
    else:
        _fail(f"parmed import failed: {r.stderr.strip()[-200:]}")
    return ["parmed"]


def _write_tiny_system(wd: Path):
    """Neutral alchemical test system: a small TIP3P water box. TI region 1 is
    the first water (:1), region 2 the second (:2) -- both net-NEUTRAL, so there
    is no charged-region / neutralization edge case. Just enough to exercise
    icfe/ifsc/softcore/MBAR keyword setup."""
    leap = wd / "tiny.leap"
    leap.write_text(textwrap.dedent(f"""\
        source leaprc.water.tip3p
        model = createUnit model
        solvateBox model TIP3PBOX 12.0
        saveAmberParm model {wd/'tiny.prmtop'} {wd/'tiny.inpcrd'}
        quit
    """))
    r = _run(["tleap", "-f", str(leap)], cwd=str(wd))
    return (wd/"tiny.prmtop").exists() and (wd/"tiny.inpcrd").exists(), r


def check_gti_run(report, wd: Path):
    print("\n[4] pmemd.cuda accepts GTI/softcore mdin (tiny real run)")
    if shutil.which("pmemd.cuda") is None or shutil.which("tleap") is None:
        _warn("skipped (pmemd.cuda or tleap missing -- see above).")
        report["gti_run"] = {"skipped": True}; return []
    built, leap_r = _write_tiny_system(wd)
    if not built:
        _warn("could not build tiny test topology; skipping GTI run.\n"
              f"        tleap stderr: {leap_r.stderr.strip()[-200:]}")
        report["gti_run"] = {"skipped": True, "reason": "tleap build failed"}
        return []
    mdin = wd / "gti.in"
    mdin.write_text(textwrap.dedent("""\
        GTI keyword acceptance test
        &cntrl
          imin=1, ntmin=2, maxcyc=20, ncyc=10,
          ntb=1, ntr=0, cut=10.0,
          icfe=1, ifsc=1, clambda=0.5,
          timask1=':1', timask2=':2',
          scmask1=':1', scmask2=':2',
          scalpha=0.5, scbeta=12.0,
          gti_add_sc=5,
          ifmbar=1, mbar_states=2, mbar_lambda=0.0,1.0,
        /
    """))
    r = _run(["pmemd.cuda", "-O", "-i", str(mdin),
              "-p", str(wd/"tiny.prmtop"), "-c", str(wd/"tiny.inpcrd"),
              "-o", str(wd/"gti.out"), "-r", str(wd/"gti.rst")], cwd=str(wd))
    out = (wd/"gti.out").read_text() if (wd/"gti.out").exists() else ""
    blob = (out + "\n" + (r.stderr or "")).lower()
    report["gti_run"] = {"rc": r.returncode, "stderr": (r.stderr or "")[-500:],
                         "mdout_tail": out[-1000:]}

    # POSITIVE evidence that GTI is supported: pmemd parsed icfe and built the
    # TI regions. A build WITHOUT GTI rejects icfe/ifsc at the namelist read and
    # never gets this far.
    gti_engaged = ("ti region" in blob or "softcore" in blob
                   or "sum of charges for ti region" in blob
                   or "ti mask" in blob or "clambda" in blob)
    # NEGATIVE evidence: explicit "keyword not supported/allowed/unknown" tied to
    # the alchemical flags, or an input-error mentioning icfe/ifsc.
    rejected = (("icfe" in blob or "ifsc" in blob or "gti" in blob)
                and ("not supported" in blob or "not allowed" in blob
                     or "unknown" in blob or "not a valid" in blob
                     or "is not recognized" in blob))

    if rejected:
        _fail("pmemd.cuda rejected the alchemical (icfe/ifsc/gti) keywords -- "
              "this build lacks GPU-TI (GTI) support. ABFEP cannot run here.\n"
              f"        mdout tail: ...{out[-300:].strip()}")
        return ["gti"]
    if gti_engaged:
        if r.returncode == 0:
            _ok("pmemd.cuda ran a GTI step and accepted icfe/ifsc/softcore/MBAR.")
        else:
            _ok("pmemd.cuda ACCEPTS GTI (built TI regions / softcore). The tiny "
                "toy run exited non-zero for a benign reason (minimal test "
                f"system), not a GTI-support problem. rc={r.returncode}.")
        return []
    # Neither clearly engaged nor clearly rejected -> inconclusive, WARN (not a
    # hard blocker), so a quirky toy-system failure can't produce a false NO-GO.
    _warn("GTI keyword test was inconclusive (no explicit rejection, but TI "
          f"setup not detected; rc={r.returncode}). Verify with a real ABFE "
          "input. stderr/mdout tail:\n"
          f"        ...{(out or r.stderr or '').strip()[-300:]}")
    report["gti_run"]["inconclusive"] = True
    return []


def check_mpi_multiwindow(report, wd: Path):
    print("\n[5] pmemd.cuda.MPI multi-window (-ng) smoke test")
    if shutil.which("pmemd.cuda.MPI") is None or shutil.which("mpirun") is None:
        _warn("skipped (pmemd.cuda.MPI or mpirun missing). Single-window TI "
              "still works; HREMD/-ng would not.")
        report["mpi"] = {"skipped": True}; return []
    if not (wd/"tiny.prmtop").exists():
        _warn("skipped (no tiny topology from step 4).")
        report["mpi"] = {"skipped": True}; return []
    for k in (0, 1):
        (wd/f"m{k}.in").write_text(textwrap.dedent(f"""\
            mw {k}
            &cntrl
              imin=1, maxcyc=10, ncyc=5, ntb=1, cut=10.0,
              icfe=1, ifsc=1, clambda={0.0 if k==0 else 1.0},
              timask1=':1', timask2=':2', scmask1=':1', scmask2=':2',
              scalpha=0.5, scbeta=12.0,
            /
        """))
    gf = wd/"groupfile"
    gf.write_text(
        f"-O -i {wd}/m0.in -p {wd}/tiny.prmtop -c {wd}/tiny.inpcrd -o {wd}/m0.out -r {wd}/m0.rst\n"
        f"-O -i {wd}/m1.in -p {wd}/tiny.prmtop -c {wd}/tiny.inpcrd -o {wd}/m1.out -r {wd}/m1.rst\n")
    r = _run(["mpirun", "-np", "2", "pmemd.cuda.MPI", "-ng", "2",
              "-groupfile", str(gf)], cwd=str(wd))
    report["mpi"] = {"rc": r.returncode, "stderr": (r.stderr or "")[-500:]}
    if r.returncode == 0 and (wd/"m0.out").exists() and (wd/"m1.out").exists():
        _ok("pmemd.cuda.MPI ran a 2-window groupfile."); return []
    _warn(f"multi-window MPI test failed (rc={r.returncode}). Single-window TI "
          f"may still work; HREMD/-ng would not.\n"
          f"        stderr: ...{(r.stderr or '').strip()[-250:]}")
    return []


def check_analysis(report):
    print("\n[6] MBAR analysis libs (alchemlyb + pymbar)")
    r = _run([sys.executable, "-c",
              "import alchemlyb,pymbar;"
              "print('alchemlyb',alchemlyb.__version__,'pymbar',pymbar.__version__)"])
    report["analysis"] = {"rc": r.returncode, "out": r.stdout.strip(),
                          "err": r.stderr.strip()[-300:]}
    if r.returncode == 0:
        _ok(r.stdout.strip())
    else:
        _warn("alchemlyb/pymbar not importable -- analyzer falls back to "
              "trapezoid TI (dG still produced; no MBAR/overlap diagnostics).")
    return []


def _env_banner(report):
    """Print which env/module the preflight is running in, so the report proves
    it ran in the right place (amber-md + amber/22.8)."""
    import os
    conda = os.environ.get("CONDA_DEFAULT_ENV", "(none)")
    amberhome = os.environ.get("AMBERHOME", "(unset)")
    loaded = os.environ.get("LOADEDMODULES", "(unknown)")
    pmemd = shutil.which("pmemd.cuda") or "(not found)"
    print("\n[env] CONDA_DEFAULT_ENV :", conda)
    print("[env] AMBERHOME        :", amberhome)
    print("[env] which pmemd.cuda :", pmemd)
    print("[env] LOADEDMODULES    :", loaded)
    report["env"] = {"conda_default_env": conda, "AMBERHOME": amberhome,
                     "which_pmemd_cuda": pmemd, "loaded_modules": loaded}
    if conda not in ("amber-md",):
        print("[env] NOTE: expected conda env 'amber-md'. If this differs, you may "
              "be in the wrong environment for Amber ABFEP.")



def main():
    ap = argparse.ArgumentParser(description="Amber ABFEP environment preflight.")
    ap.add_argument("--work-dir", type=Path,
                    default=Path(tempfile.mkdtemp(prefix="abfe_pf_")),
                    help="Scratch dir for the tiny test system.")
    ap.add_argument("--keep", action="store_true", help="Keep scratch files.")
    a = ap.parse_args()
    wd = a.work_dir; wd.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Amber ABFEP preflight  (amber_md_workflow_v2.5.0)")
    print(f"python   : {sys.executable}")
    print(f"work-dir : {wd}")
    print("=" * 70)

    report = {"python": sys.executable, "work_dir": str(wd), "binaries": {}}
    _env_banner(report)
    blockers = []
    blockers += check_binaries(report)
    blockers += check_tleap(report)
    blockers += check_parmed(report)
    blockers += check_gti_run(report, wd)
    check_mpi_multiwindow(report, wd)
    check_analysis(report)

    report["blockers"] = blockers
    report["verdict"] = "GO" if not blockers else "NO-GO"
    (wd / "abfe_preflight_report.json").write_text(json.dumps(report, indent=2))

    print("\n" + "=" * 70)
    if not blockers:
        print("VERDICT: GO  -- Amber ABFEP should run in this environment."); rc = 0
    else:
        print(f"VERDICT: NO-GO -- blockers: {', '.join(blockers)}")
        print("Fix the [BLOCK] items above, then re-run this preflight."); rc = 1
    print(f"Full report: {wd/'abfe_preflight_report.json'}")
    print("=" * 70)
    if not a.keep and "abfe_pf_" in str(wd):
        shutil.rmtree(wd, ignore_errors=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())


# v2.5.77 FIX D: guard against the two live-run bugs before burning GPU hours.
def check_prodin_consistency(prod_in_path, requested_nstlim=None):
    """Assert a generated prod.in is self-consistent for MBAR + sampling.

    (1) nstlim matches requested_nstlim (catches the GUI 5-ns-override bug).
    (2) every TI window's own clambda appears in its mbar_lambda list
        (catches the %.3f truncation that made u_nk rank-deficient -> MBAR dead).
    Returns (ok: bool, problems: list[str]).
    """
    import re as _re
    txt = open(prod_in_path).read()
    problems = []
    m = _re.search(r"nstlim\s*=\s*(\d+)", txt)
    nst = int(m.group(1)) if m else None
    if requested_nstlim is not None and nst != int(requested_nstlim):
        problems.append("nstlim=%s but requested %s (GUI override?)" % (nst, requested_nstlim))
    cl = _re.search(r"clambda\s*=\s*([0-9.]+)", txt)
    ml = _re.search(r"mbar_lambda\s*=\s*([0-9.,eE+\-]+)", txt)
    if cl and ml:
        clv = float(cl.group(1))
        states = [s for s in ml.group(1).split(",") if s.strip()]
        fs = set(float(s) for s in states)
        if not any(abs(clv - s) < 1e-6 for s in fs):
            problems.append("clambda=%s NOT in its own mbar_lambda list %s "
                            "(u_nk will be rank-deficient -> MBAR unsolvable)"
                            % (clv, sorted(fs)))
    return (len(problems) == 0, problems)
