"""amber_md.mdin_validator  (v2.5.31)

Importable core for the GPU-free pre-submit mdin checks. Both the CLI wrapper
(tools/validate_mdin.py) and the driver (amber_md.fep_driver) call into here, so
there is a single source of truth for the rules.

It catches the *setup-time parser* failures that have repeatedly aborted runs only
AFTER reaching a GPU node -- one defect at a time, because pmemd aborts on the
first error it meets:

  * inline "&wt type=..." weight-change card  -> "Invalid TYPE flag in line"   (v2.5.29)
  * "*" residue wildcard in a restraintmask    -> "unknown symbol:*"            (v2.5.28)
  * ntr=1 stage whose run command lacks -ref    -> "Unit 10 Error on OPEN: refc" (v2.5.27)

plus heat-is-NVT, icfe=1 TI-mask presence, and DISANG-file-exists checks.

Public API:
    validate_leg(leg_dir)      -> (errors:list[Issue], warnings:list[Issue])
    validate_paths(paths)      -> (errors, warnings, n_windows)
    format_issues(errs, warns) -> str

Only the standard library is used.
"""
from __future__ import annotations
import re
from pathlib import Path

RESTRAINED_STAGES = {"heat", "dens", "eq"}
ALL_STAGES = ["min", "heat", "dens", "eq", "prod"]
REF_SOURCE = {"heat": "min.rst", "dens": "heat.rst (or min.rst)", "eq": "dens.rst"}


class Issue:
    __slots__ = ("level", "where", "msg")

    def __init__(self, level: str, where: str, msg: str):
        self.level = level
        self.where = where
        self.msg = msg

    def __str__(self) -> str:
        tag = "ERROR" if self.level == "ERROR" else "warn "
        return f"  [{tag}] {self.where}: {self.msg}"


def _read(p: Path) -> str:
    return p.read_text(errors="replace")


def _cntrl_flag(text: str, name: str):
    m = re.search(rf"\b{name}\s*=\s*([-\w.'\"]+)", text)
    return m.group(1).strip().strip("'\"") if m else None


_WILDCARD_RES = re.compile(r":[^&|@]*\*")



WT_CANON_TYPES = {
    "BOND","ANGLE","TORSION","IMPROP","RSTAR","INTERN","ALL","ELEC","NB","HB",
    "ATTRACT","REPULSE","VDW","SOFTREPULSION","EPB","RESTRAINT","REST",
    "RESTRAINT_WT","TGT","TGTMDFRC","TEMP0","TEMP0LES","TAUTP","CUT","NSTEP0",
    "NSTEP1","STPMLT","NPRT","DISAVE","ANGAVE","TORAVE","DISAVI","ANGAVI",
    "TORAVI","DIPMASS","DIPSCALE","DUMPAVE","DUMPFREQ","END",
}
WT_TYPO_HINT = {"TEMP_0":"TEMP0","TEMP_0LES":"TEMP0LES","RESTRAIN":"RESTRAINT","RESTRAINT_W":"RESTRAINT_WT"}

def check_wt_type_keywords(text: str, where: str, issues: list) -> None:
    for m in re.finditer(r"type\s*=\s*['\"]?([A-Za-z0-9_]+)['\"]?", text):
        tok=m.group(1).upper()
        if tok in WT_CANON_TYPES: continue
        hint=WT_TYPO_HINT.get(tok)
        if hint:
            issues.append(Issue("ERROR",where,f"&wt type='{m.group(1)}' is not a valid Amber keyword -> 'Invalid TYPE flag' aborts pmemd at heat. Use type='{hint}'."))
        else:
            issues.append(Issue("ERROR",where,f"&wt type='{m.group(1)}' is not a recognized Amber weight-change keyword -> pmemd aborts with 'Invalid TYPE flag'. Valid: TEMP0, RESTRAINT, DUMPFREQ, END."))


def check_wt_cards(text: str, where: str, issues: list) -> None:
    """Validate the &wt weight-change section.

    pmemd's nmropt=1 reader requires a DATA-bearing &wt card (TEMP_0, DUMPFREQ,
    REST, ...) to have its '&wt' opener alone on its line; an inline data card is
    mis-tokenized -> 'Invalid TYPE flag' (the v2.5.29 heat regression).

    EXCEPTION (v2.5.31a): the bare terminator written inline as "&wt type='END' /"
    is TOLERATED by pmemd and has run in the decharge/vdw/restraint dens/eq/prod
    tails for years (it carries no istep/value data to mis-read). We must NOT flag
    it, or the gate would block legs that demonstrably work. Only a DATA card --
    i.e. an inline &wt whose type is not END, or which carries istep/value fields
    -- is fatal.
    """
    if _cntrl_flag(text, "nmropt") != "1":
        if "&wt" in text:
            issues.append(Issue("WARN", where, "&wt block present but nmropt != 1"))
        return
    for ln in text.splitlines():
        s = ln.strip()
        if not s.startswith("&wt") or s == "&wt":
            continue
        # inline &wt ... -> is it the harmless END terminator, or a data card?
        is_end_only = bool(re.match(r"&wt\s+type\s*=\s*'END'\s*/?\s*$", s))
        if is_end_only:
            continue  # proven-safe inline terminator; do not flag
        issues.append(Issue(
            "ERROR", where,
            f"inline &wt DATA card -> pmemd 'Invalid TYPE flag'. A data-bearing "
            f"&wt card (e.g. TEMP_0) needs '&wt' alone on its line, then "
            f"'type = ...' on the next. Offending line: {s!r}"))
    if "&wt" in text and "type = 'END'" not in text and "type='END'" not in text:
        issues.append(Issue("ERROR", where,
                            "nmropt=1 with &wt cards but no &wt type='END' terminator"))


def check_masks(text: str, where: str, issues: list) -> None:
    for key in ("restraintmask", "bellymask"):
        m = re.search(rf"{key}\s*=\s*['\"]([^'\"]*)['\"]", text)
        if m:
            val = m.group(1)
            if "*" in val and _WILDCARD_RES.search(val):
                issues.append(Issue(
                    "ERROR", where,
                    f"{key}={val!r} contains a '*' residue wildcard; pmemd's "
                    f"group-input parser rejects it ('unknown symbol:*'). Use "
                    f"explicit residue names (e.g. WAT,HOH,Na+,Cl-,K+)."))
    if _cntrl_flag(text, "icfe") == "1":
        # timask1 must be present. timask2 is REQUIRED only for a two-region
        # transformation (RBFE/relative: ligand A <-> ligand B). For ABFE
        # (absolute decoupling of a single molecule) timask2 is LEGITIMATELY
        # empty -- e.g. the driver logs "timask1=':LIG' timask2=''". Flagging an
        # empty timask2 there is a false positive that blocks valid ABFE runs.
        m1 = re.search(r"timask1\s*=\s*['\"]([^'\"]*)['\"]", text)
        if not m1 or not m1.group(1).strip():
            issues.append(Issue("ERROR", where,
                                "icfe=1 but timask1 is missing/empty"))
        # Heuristic for 'this is RBFE': a non-empty scmask2/crgmask2 or two
        # distinct residues in the masks would imply a second region. We only
        # require timask2 when scmask2 is present and non-empty.
        m_sc2 = re.search(r"scmask2\s*=\s*['\"]([^'\"]*)['\"]", text)
        if m_sc2 and m_sc2.group(1).strip():
            m2 = re.search(r"timask2\s*=\s*['\"]([^'\"]*)['\"]", text)
            if not m2 or not m2.group(1).strip():
                issues.append(Issue(
                    "ERROR", where,
                    "icfe=1 with a non-empty scmask2 (two-region/RBFE) but "
                    "timask2 is missing/empty"))


def check_restraint_internal(text: str, where: str, issues: list) -> bool:
    """Returns True if this stage is restrained (ntr=1) and thus needs -ref."""
    if _cntrl_flag(text, "ntr") == "1":
        rm = re.search(r"restraintmask\s*=\s*['\"]([^'\"]*)['\"]", text)
        if not rm or not rm.group(1).strip():
            issues.append(Issue("ERROR", where, "ntr=1 but restraintmask is missing/empty"))
        return True
    return False


def check_heat_is_nvt(text: str, where: str, issues: list) -> None:
    if _cntrl_flag(text, "ntp") not in (None, "0"):
        issues.append(Issue("ERROR", where, f"heat must be NVT but ntp={_cntrl_flag(text,'ntp')}"))
    if _cntrl_flag(text, "ntb") not in (None, "1"):
        issues.append(Issue("WARN", where, f"heat expected ntb=1, got ntb={_cntrl_flag(text,'ntb')}"))


def check_disang(text: str, wd: Path, where: str, issues: list) -> None:
    m = re.search(r"^\s*DISANG\s*=\s*(\S+)", text, re.M)
    if m:
        fn = m.group(1).strip()
        if not (wd / fn).exists():
            issues.append(Issue("ERROR", where, f"DISANG={fn} but {fn} not found in {wd.name}/"))


def check_run_script(leg: Path, needref_stages: set, issues: list,
                     require_script: bool = True) -> None:
    """Verify the leg's run script passes -ref for every restrained stage.

    The driver writes ONE script per leg, named run_<leg>.lsf, in the leg dir
    (the parent of the lambda_* windows). We must NOT treat analyze_*.lsf or
    cycle_close_*.lsf as the MD run script. Called ONCE per leg.

    require_script: when False (validator runs before the script is written, e.g.
    the driver's pre-submit content gate), a missing script is silently ignored
    instead of warned -- the post-submit pass does the real -ref check.
    """
    cand = sorted(p for p in leg.glob("run_*.lsf"))
    if not cand:
        # fall back to any .lsf that is not an analyze/cycle-closer script
        cand = sorted(p for p in leg.glob("*.lsf")
                      if not p.name.startswith(("analyze_", "cycle_close")))
    if not cand:
        cand = sorted(leg.glob("run*.sh")) + sorted(leg.glob("submit*.sh"))
    if not cand:
        if require_script:
            issues.append(Issue("WARN", leg.name,
                                "no run_<leg>.lsf found to verify -ref"))
        return
    rs = cand[0]
    txt = _read(rs)
    for stage in sorted(needref_stages):
        # match the stage's pmemd line in either format:
        #   run_stage heat pmemd.cuda -O -i heat.in ... -ref min.rst
        #   pmemd.cuda -O -i heat.in ... -ref ...
        m = re.search(rf"-i\s+{stage}\.in\b[^\n]*", txt)
        if not m:
            issues.append(Issue("WARN", rs.name, f"no run line found for {stage}.in"))
            continue
        if "-ref" not in m.group(0):
            issues.append(Issue(
                "ERROR", rs.name,
                f"{stage} runs ntr=1 (restrained) but its pmemd command has no "
                f"-ref -> 'Unit 10 Error on OPEN: refc'. Add "
                f"'-ref {REF_SOURCE.get(stage,'<start coords>')}'."))
def _validate_window(wd: Path, issues: list) -> set:
    """Run all window-level mdin checks. Returns the set of restrained stages
    (those needing -ref) so the caller can do ONE leg-level run-script check.
    Does NOT touch the run script itself (that is a per-leg concern)."""
    where_pref = f"{wd.parent.name}/{wd.name}"
    present = [s for s in ALL_STAGES if (wd / f"{s}.in").exists()]
    if not present:
        return set()
    needref = set()
    for stage in present:
        text = _read(wd / f"{stage}.in")
        where = f"{where_pref}/{stage}.in"
        check_wt_cards(text, where, issues)
        check_wt_type_keywords(text, where, issues)
        check_masks(text, where, issues)
        check_disang(text, wd, where, issues)
        if stage == "heat":
            check_heat_is_nvt(text, where, issues)
        if stage in RESTRAINED_STAGES and check_restraint_internal(text, where, issues):
            needref.add(stage)
    return needref


def _find_windows(root: Path):
    if (root / "min.in").exists() or (root / "heat.in").exists():
        yield root
        return
    for lam in sorted(root.glob("lambda_*")):
        if lam.is_dir():
            yield lam
    for sub in sorted(root.glob("*")):
        if sub.is_dir() and not sub.name.startswith("lambda_"):
            for lam in sorted(sub.glob("lambda_*")):
                if lam.is_dir():
                    yield lam


def validate_leg(leg_dir, require_script: bool = True):
    """Validate one leg directory (lambda_* windows + run_<leg>.lsf).
    Returns (errors, warnings).

    require_script=False suppresses the 'no run_<leg>.lsf' warning -- use it for the
    driver's pre-submit content gate (the script is written later by submit_leg);
    the driver then calls validate_leg again with require_script=True after submit
    to actually verify -ref.
    """
    issues: list = []
    leg = Path(leg_dir)
    needref: set = set()
    for wd in _find_windows(leg):
        needref |= _validate_window(wd, issues)
    if needref:
        check_run_script(leg, needref, issues, require_script=require_script)
    errs = [i for i in issues if i.level == "ERROR"]
    warns = [i for i in issues if i.level == "WARN"]
    return errs, warns

def validate_paths(paths):
    """Validate one or more leg/run directories. Returns (errors, warnings, n_windows)."""
    issues: list = []
    seen = set()
    n = 0
    for p in paths:
        root = Path(p).expanduser()
        if not root.exists():
            issues.append(Issue("ERROR", str(root), "path does not exist"))
            continue
        # group windows by their parent leg so the run-script check runs once per leg
        legs: dict = {}
        for wd in _find_windows(root):
            key = str(wd.resolve())
            if key in seen:
                continue
            seen.add(key)
            n += 1
            needref = _validate_window(wd, issues)
            legs.setdefault(str(wd.parent.resolve()), (wd.parent, set()))
            legs[str(wd.parent.resolve())][1].update(needref)
        for _k, (legdir, needref) in legs.items():
            if needref:
                check_run_script(legdir, needref, issues, require_script=True)
    errs = [i for i in issues if i.level == "ERROR"]
    warns = [i for i in issues if i.level == "WARN"]
    return errs, warns, n


def format_issues(errs, warns) -> str:
    """Render Issue lists to a printable block (errors first, then warnings).

    v2.5.31f: restored. The body had been orphaned as dead code after the
    `return` in validate_paths(), leaving format_issues undefined -> fep_driver
    and tools/validate_mdin.py failed to import it and the mdin gate was silently
    skipped ("could not run ... cannot import name 'format_issues'; continuing").
    """
    lines = []
    for i in list(errs) + list(warns):
        lines.append(str(i))
    return "\n".join(lines)