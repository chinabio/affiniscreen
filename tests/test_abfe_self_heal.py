"""Tests for amber_md/abfe_self_heal.py against the REAL detonated eq.out.

Uses the same fixtures as the estimator regression test:
    fixtures/eq_complex_vdw_lambda_0.500.out  -- detonated (14,963 K, **** MBAR)

Run:  pytest -q test_abfe_self_heal.py
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import os
import shutil
import textwrap
import pytest

import abfe_self_heal as sh

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
EQ_BLOWN = os.path.join(FIX, "eq_complex_vdw_lambda_0.500.out")


_MIN_MDIN = """\
FEP prod, lambda=0.5
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=1000000, dt=0.002,
  ntc=2, ntf=1, ntt=3, gamma_ln=2.0,
  tempi=100.0, temp0=298.0,
  ntp=1, pres0=1.0, taup=2.0,
  icfe=1, ifsc=1, clambda=0.5,
  scalpha=0.5, scbeta=12.0,
  nscm=0,
/
"""


def _make_window(tmp_path, with_blown_eq=True, prod_complete=False,
                 err_text=None, lsf_out=None):
    wd = tmp_path / "lambda_0.500"
    wd.mkdir()
    for f in ("eq.in", "prod.in", "dens.in"):
        (wd / f).write_text(_MIN_MDIN)
    if with_blown_eq:
        shutil.copy(EQ_BLOWN, wd / "eq.out")
    if prod_complete:
        (wd / "prod.out").write_text(
            "NSTEP = 1000000\n A V E R A G E S   O V E R 1000000 STEPS\n"
            "Final Performance Info\n")
    if err_text is not None:
        (wd / "fep_complex_vdw.7.123.err").write_text(err_text)
    if lsf_out is not None:
        (wd / "fep_complex_vdw.7.123.out").write_text(lsf_out)
    return wd
    return wd


# --------------------------------------------------------------------------- #
# completion detection                                                         #
# --------------------------------------------------------------------------- #
def test_prod_complete_true(tmp_path):
    wd = _make_window(tmp_path, prod_complete=True)
    assert sh.prod_is_complete(wd / "prod.out") is True


def test_prod_complete_false_when_missing(tmp_path):
    wd = _make_window(tmp_path, prod_complete=False)
    assert sh.prod_is_complete(wd / "prod.out") is False


def test_prod_complete_false_when_truncated(tmp_path):
    wd = _make_window(tmp_path)
    (wd / "prod.out").write_text("NSTEP = 5000\n (run was killed mid-way)\n")
    assert sh.prod_is_complete(wd / "prod.out") is False


# --------------------------------------------------------------------------- #
# classification on the REAL detonated eq.out                                  #
# --------------------------------------------------------------------------- #
def test_classify_real_blowup(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=True)
    diag = sh.classify_failure(wd)
    assert diag.failure_class == "blowup_temperature"
    assert diag.evidence  # carries the matched evidence line


def test_classify_box_drift(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=False)
    (wd / "prod.out").write_text(
        "Started prod...\n"
        "Periodic box dimensions have changed too much from their initial values\n")
    diag = sh.classify_failure(wd)
    assert diag.failure_class == "box_drift"


def test_classify_oom_from_err(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=False,
                      err_text="TERM_MEMLIMIT: job killed after reaching mem limit")
    diag = sh.classify_failure(wd)
    assert diag.failure_class == "oom_or_killed"


# Real LSF job-summary tail from the gpu-node-02 exit-141 failure (complex_vdw).
_REAL_EXIT141_LSF_OUT = (
    "[run_stage] lambda=0.300 stage=eq starting: pmemd.cuda -O -i eq.in ...\n"
    "Sender: LSF System <lsfadmin@gpu-node-02>\n"
    "Subject: Job 4694658[7]: <fep_complex_vdw[1-23]> ... Exited\n"
    "Started at Wed Jun 10 20:56:35 2026\n"
    "Terminated at Wed Jun 10 21:24:11 2026\n"
    "Exited with exit code 141.\n"
)


def test_classify_external_kill_exit141_from_real_lsf_out(tmp_path):
    """The real failure: pmemd ran eq, then LSF terminated the job (exit 141).
    Only IEEE_UNDERFLOW in .err; the verdict is in the LSF .out summary."""
    wd = _make_window(tmp_path, with_blown_eq=False,
                      err_text=("Note: The following floating-point exceptions "
                                "are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL\n"),
                      lsf_out=_REAL_EXIT141_LSF_OUT)
    diag = sh.classify_failure(wd)
    assert diag.failure_class == "external_kill", diag
    assert "141" in diag.evidence or "Terminated" in diag.evidence


def test_classify_external_kill_term_runlimit(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=False,
                      lsf_out="TERM_RUNLIMIT: job killed after reaching LSF run "
                              "time limit.\nExited with exit code 140.\n")
    assert sh.classify_failure(wd).failure_class == "external_kill"


def test_policy_external_kill_is_clean_resubmit():
    p1 = sh.remediation_plan("external_kill", 1)
    assert "__resubmit__" in p1 and p1["__resubmit__"]["edit_mdin"] is False
    # no mdin file edits on a clean resubmit
    assert not any(k.endswith(".in") for k in p1)
    # escalates to host-avoidance on attempt 2
    p2 = sh.remediation_plan("external_kill", 2)
    assert p2.get("__resources__", {}).get("avoid_last_host") is True

def test_classify_missing_no_error(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=False)
    diag = sh.classify_failure(wd)
    assert diag.failure_class == "missing_no_error"


# --------------------------------------------------------------------------- #
# remediation policy escalates                                                 #
# --------------------------------------------------------------------------- #
def test_policy_blowup_shrinks_dt_and_escalates():
    p1 = sh.remediation_plan("blowup_temperature", 1)
    p2 = sh.remediation_plan("blowup_temperature", 2)
    p3 = sh.remediation_plan("blowup_temperature", 3)
    assert p1["eq.in"]["dt"] == 0.001
    assert p2["eq.in"]["dt"] == 0.0005
    assert p2["eq.in"]["gamma_ln"] >= p1["eq.in"]["gamma_ln"]
    # last resort softens the soft-core core
    assert p3["eq.in"]["scbeta"] > 12.0


def test_policy_softcore_com_disables_nscm():
    p = sh.remediation_plan("softcore_com", 1)
    assert p["eq.in"]["nscm"] == 0 and p["prod.in"]["nscm"] == 0


def test_policy_oom_requests_resources():
    p = sh.remediation_plan("oom_or_killed", 1)
    assert "__resources__" in p and p["__resources__"]["fep_mem_mb"] > 8192


# --------------------------------------------------------------------------- #
# mdin editing                                                                 #
# --------------------------------------------------------------------------- #
def test_apply_mdin_edits_changes_existing_and_inserts_missing(tmp_path):
    f = tmp_path / "eq.in"
    f.write_text(_MIN_MDIN)
    changes = sh.apply_mdin_edits(f, {"dt": 0.0005, "gamma_ln": 10.0,
                                      "barostat": 1})
    txt = f.read_text()
    assert "dt=0.0005," in txt.replace(" ", "")
    assert "gamma_ln=10," in txt.replace(" ", "")
    assert "barostat=1," in txt.replace(" ", "")          # inserted (was absent)
    assert (tmp_path / "eq.in.orig").exists()             # backup kept
    names = {c[0] for c in changes}
    assert {"dt", "gamma_ln", "barostat"} <= names


def test_apply_mdin_edits_verbatim_shell_ref(tmp_path):
    f = tmp_path / "eq.in"
    f.write_text(_MIN_MDIN)
    sh.apply_mdin_edits(f, {"tempi": "${TEMP0}"})
    assert "tempi=${TEMP0}," in f.read_text().replace(" ", "")


# --------------------------------------------------------------------------- #
# end-to-end heal loop with a stubbed rerun                                    #
# --------------------------------------------------------------------------- #
def test_heal_window_converges_after_remediation(tmp_path):
    """First the window has a detonated eq and no prod.out. The stubbed rerun
    'succeeds' (writes a complete prod.out) only AFTER dt has been reduced,
    emulating the real fix (gentle integrator -> stable run)."""
    wd = _make_window(tmp_path, with_blown_eq=True)
    calls = {"n": 0}

    def rerun(window_dir, attempt):
        calls["n"] += 1
        eq_txt = (window_dir / "eq.in").read_text().replace(" ", "")
        # only "stabilises" once dt was shrunk to <= 0.001
        if "dt=0.001," in eq_txt or "dt=0.0005," in eq_txt:
            (window_dir / "prod.out").write_text(
                "A V E R A G E S   O V E R 1000000 STEPS\nFinal Performance Info\n")
            return 0
        return 1

    rep = sh.heal_window(wd, rerun, max_attempts=3)
    assert rep.completed is True
    assert rep.attempts == 1                  # first remediation already sets dt=0.001
    assert rep.history[0]["failure_class"] == "blowup_temperature"
    assert rep.history[0]["complete_after"] is True


def test_heal_window_gives_up_and_reports(tmp_path):
    """If reruns never succeed, the loop exhausts attempts and reports cleanly."""
    wd = _make_window(tmp_path, with_blown_eq=True)

    def rerun_always_fail(window_dir, attempt):
        return 1

    rep = sh.heal_window(wd, rerun_always_fail, max_attempts=3)
    assert rep.completed is False
    assert rep.attempts == 3
    assert len(rep.history) == 3


def test_heal_window_noop_when_already_complete(tmp_path):
    wd = _make_window(tmp_path, with_blown_eq=True, prod_complete=True)

    def rerun_should_not_be_called(window_dir, attempt):
        raise AssertionError("rerun must not be called when already complete")

    rep = sh.heal_window(wd, rerun_should_not_be_called, max_attempts=3)
    assert rep.completed is True and rep.attempts == 0


def test_heal_leg_covers_all_windows(tmp_path):
    """heal_leg applies uniformly across a leg's windows (any of the 4 legs)."""
    leg = tmp_path / "complex_vdw"
    leg.mkdir()
    for lam in (0.0, 0.5, 1.0):
        wd = leg / f"lambda_{lam:.3f}"
        wd.mkdir()
        for f in ("eq.in", "prod.in", "dens.in"):
            (wd / f).write_text(_MIN_MDIN)
        # make 0.5 already-complete, others missing
        if lam == 0.5:
            (wd / "prod.out").write_text("Final Performance Info\n")

    def rerun(window_dir, attempt):
        (window_dir / "prod.out").write_text("Final Performance Info\n")
        return 0

    reps = sh.heal_leg(leg, (0.0, 0.5, 1.0), rerun, max_attempts=2)
    assert set(reps) == {"0.000", "0.500", "1.000"}
    assert all(r.completed for r in reps.values())
    assert reps["0.500"].attempts == 0       # was already complete


# --------------------------------------------------------------------------- #
# v2.5.8: true host-avoidance
# --------------------------------------------------------------------------- #
import abfe_self_heal_cli as cli


def test_failed_hosts_parsed_from_lsf_out(tmp_path):
    leg = tmp_path / "complex_decharge"
    wd = leg / "lambda_0.500"
    wd.mkdir(parents=True)
    (leg / "fep_complex_decharge.8.4694666.out").write_text(
        "Job was executed on host(s) <gpu-node-02>, in queue <gpu> ...\n"
        "Exited with exit code 141.\n")
    (wd / "heal_lambda_0.500.5001.out").write_text(
        "Job was executed on host(s) <1*gpu-node-01> ...\nExited with exit code 141.\n")
    (wd / "eq.out").write_text("NSTEP = 1 executed on host(s) <bogus>\n")
    hosts = cli._failed_hosts(wd)
    assert hosts == {"gpu-node-02", "gpu-node-01"}, hosts


def test_select_resource_expression():
    assert cli._select_resource(set()) == ""
    line = cli._select_resource({"gpu-node-02", "gpu-node-01"})
    assert line == ('#BSUB -R "select[hname!=\'gpu-node-01\' && '
                    'hname!=\'gpu-node-02\']"')


def test_bsub_rerun_script_excludes_failed_and_blocklist_hosts(tmp_path, monkeypatch):
    leg = tmp_path / "complex_decharge"
    wd = leg / "lambda_0.500"
    wd.mkdir(parents=True)
    (wd / "eq.rst").write_text("x")
    (leg / "fep_complex_decharge.8.4694666.out").write_text(
        "Job was executed on host(s) <gpu-node-02> ...\nExited with exit code 141.\n")

    captured = {}
    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        class R: returncode = 0
        return R()
    monkeypatch.setattr(cli.subprocess, "run", fake_run)

    rerun = cli.make_bsub_rerun("proj", "gpu", "24:00", 8192,
                                avoid_hosts=["gpu-node-01"])
    rc = rerun(wd, attempt=2)
    assert rc == 0
    script = (wd / "run_heal.lsf").read_text()
    assert "select[" in script
    assert "hname!='gpu-node-01'" in script
    assert "hname!='gpu-node-02'" in script
    assert 'cd "' in script and "module load amber/22" in script
    assert captured["cmd"] == ["bsub", "-K"]


def test_no_host_avoidance_still_honors_blocklist(tmp_path, monkeypatch):
    leg = tmp_path / "complex_decharge"
    wd = leg / "lambda_0.500"
    wd.mkdir(parents=True)
    (wd / "eq.rst").write_text("x")
    (leg / "fep_complex_decharge.8.4694666.out").write_text(
        "Job was executed on host(s) <gpu-node-02> ...\n")
    monkeypatch.setattr(cli.subprocess, "run",
                        lambda *a, **k: type("R", (), {"returncode": 0})())
    rerun = cli.make_bsub_rerun("proj", "gpu", "24:00", 8192,
                                avoid_hosts=["gpu-node-01"],
                                host_avoidance=False)
    rerun(wd, attempt=1)
    script = (wd / "run_heal.lsf").read_text()
    assert "hname!='gpu-node-01'" in script
    assert "hname!='gpu-node-02'" not in script