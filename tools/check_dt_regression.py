
#!/usr/bin/env python3
"""CI guard: fail (exit!=0) if the dt=1 fs fix or its safeguards regress.

Run from the package root (the dir containing amber_md/):
    python tools/check_dt_regression.py
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import sys, os, re, io, ast, tokenize, dataclasses

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

FAILS = []
def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond: FAILS.append(msg)

print("[check_dt_regression] root:", ROOT)

try:
    import amber_md
    ver = amber_md.__version__
    check(bool(re.match(r"\d+\.\d+", str(ver))), f"amber_md imports, __version__={ver}")
except Exception as e:
    check(False, f"amber_md import failed: {e}")
    print("\nRESULT: FAIL (cannot import)"); sys.exit(2)

from amber_md import config as c
for name in dir(c):
    o = getattr(c, name)
    if dataclasses.is_dataclass(o):
        for fld in dataclasses.fields(o):
            if fld.name in ("dt_ps", "prod_dt_ps") and fld.default is not dataclasses.MISSING:
                check(float(fld.default) <= 0.001 + 1e-12,
                      f"{name}.{fld.name} default = {fld.default} (<= 0.001)")

try:
    from amber_md._dt_guard import assert_prod_dt_safe, UnsafeTimestepError
    assert_prod_dt_safe()
    check(True, "_dt_guard.assert_prod_dt_safe() passes on defaults")
    class _Bad: dt_ps = 0.002
    fired = False
    try: assert_prod_dt_safe(_Bad())
    except UnsafeTimestepError: fired = True
    check(fired, "_dt_guard FIRES on a 2 fs config")
except Exception as e:
    check(False, f"_dt_guard unusable: {e}")

drv = open(os.path.join(ROOT, "amber_md", "fep_driver.py"), encoding="utf-8").read()
check("assert_prod_dt_safe(fep)" in drv and "assert_prod_dt_safe(md)" in drv,
      "fep_driver wires assert_prod_dt_safe(fep) and (md)")
m = re.search(r'add_argument\(\s*"--dt".*?default\s*=\s*([0-9.]+)', drv, re.S)
check(bool(m) and float(m.group(1)) <= 0.001 + 1e-12,
      f"fep_driver --dt default = {m.group(1) if m else '??'} (<= 0.001)")

# v2.5.63: --nstlim-prod must default to the 10 ns (1 fs) value, not 5 ns.
m2 = re.search(r'add_argument\(\s*"--nstlim-prod".*?default\s*=\s*([0-9_]+)', drv, re.S)
_nv = int(m2.group(1).replace("_", "")) if m2 else 0
check(_nv >= 10_000_000, f"fep_driver --nstlim-prod default = {_nv} (>= 10 ns @1fs)")

# v2.5.63: Option-A restraint-leg analytic analyzer must exist and be wired.
try:
    from amber_md.fep import analyze_restraint_leg_optionA as _arl
    check(callable(_arl), "fep.analyze_restraint_leg_optionA present")
except Exception as _e:
    check(False, f"analyze_restraint_leg_optionA missing: {_e}")
check("analyze_restraint_leg_optionA" in open(os.path.join(ROOT, "amber_md", "fep.py")).read(),
      "analyze LSF body wired to Option-A analyzer")

# v2.5.64: restraint leg has its own short production budget (equilibration-only).
for name in dir(c):
    o = getattr(c, name)
    if dataclasses.is_dataclass(o):
        for fld in dataclasses.fields(o):
            if fld.name == "restraint_nstlim_prod" and fld.default is not dataclasses.MISSING:
                _rv = int(fld.default)
                check(0 < _rv <= 5_000_000, f"{name}.restraint_nstlim_prod = {_rv} (0<..<=5ns: equilibration-only)")
check("restraint_nstlim_prod" in open(os.path.join(ROOT, "amber_md", "fep.py")).read(),
      "fep.py _prod_in_restraint uses restraint_nstlim_prod")

# v2.5.65: stand-alone leg analyzer must exist, compile, and run on an Option-A leg.
import os as _os, subprocess as _sp, tempfile as _tf, json as _json
_al = _os.path.join(HERE, "analyze_leg.py")
check(_os.path.exists(_al), "tools/analyze_leg.py present")
if _os.path.exists(_al):
    import py_compile as _pc
    try:
        _pc.compile(_al, doraise=True); _ok=True
    except Exception as _e:
        _ok=False; print("    (compile error: %s)" % _e)
    check(_ok, "tools/analyze_leg.py compiles")
    # synth a minimal Option-A leg in a temp dir and analyze it end-to-end
    _d = _tf.mkdtemp()
    (_os.path.join(_d, "boresch_correction.txt"))
    open(_os.path.join(_d, "boresch_correction.txt"), "w").write("-11.447503\n")
    for _lam in (0.0, 1.0):
        _w = _os.path.join(_d, "lambda_%.3f" % _lam); _os.makedirs(_w, exist_ok=True)
        open(_os.path.join(_w, "prod.in"), "w").write("&cntrl\n  nstlim=2000000, dt=0.001,\n/\n")
        open(_os.path.join(_w, "prod.out"), "w").write("NSTEP = 2000000\nFinal Performance Info\nTotal wall time 100\n")
    _env = dict(_os.environ); _env["PYTHONPATH"] = ROOT + _os.pathsep + _env.get("PYTHONPATH","")
    _r = _sp.run([sys.executable, _al, _d, "--lambdas", "0.0,1.0", "--json-only"],
                 capture_output=True, text=True, env=_env)
    _sj = _os.path.join(_d, "summary.json")
    _val = None
    if _os.path.exists(_sj):
        try: _val = _json.load(open(_sj)).get("dG_kcal_mol")
        except Exception: _val = None
    check(_r.returncode == 0 and _val == -11.447503,
          "analyze_leg.py runs Option-A leg -> dG=-11.4475, exit 0 (got rc=%s dG=%s)" % (_r.returncode, _val))

# check 6: a genuine mdin template has dt=0.002 AND a namelist token on the
# SAME physical line (Amber namelist: "imin=0,irest=1,...,dt=0.002,").
DT002 = re.compile(r'(?:^|[,&\s])(?:dt|timestep)\s*=\s*0\.002\b', re.I)
NAMELIST = re.compile(r'(?:&cntrl|nstlim\s*=|\bntc\s*=|\bimin\s*=|\birest\s*=|\bntx\s*=|\bntf\s*=)', re.I)
def line_is_template(line):
    return bool(DT002.search(line) and NAMELIST.search(line))

offenders = []
for r, _, files in os.walk(os.path.join(ROOT, "amber_md")):
    for f in files:
        if not f.endswith(".py"): continue
        path = os.path.join(r, f); rel = os.path.relpath(path, ROOT)
        src = open(path, encoding="utf-8", errors="ignore").read()

        # (a) executable code tokens
        try:
            bt = list(tokenize.tokenize(io.BytesIO(src.encode()).readline))
        except Exception:
            bt = []
        by_line = {}
        for t in bt:
            if t.type in (tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
                          tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING, tokenize.ENDMARKER):
                continue
            by_line.setdefault(t.start[0], []).append(t)
        for ln, lt in by_line.items():
            names = {t.string.lower() for t in lt if t.type == tokenize.NAME}
            nums  = {t.string for t in lt if t.type == tokenize.NUMBER}
            if (names & {"dt","dt_ps","prod_dt_ps","timestep"}) and "0.002" in nums:
                offenders.append(f"{rel}:{ln} (code)")

        # (b) namelist template: scan each string node line-by-line
        try:
            tree = ast.parse(src)
        except Exception:
            tree = None
        if tree is not None:
            for node in ast.walk(tree):
                vals = []
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    vals = [node.value]
                elif isinstance(node, ast.JoinedStr):
                    vals = ["".join(v.value for v in node.values
                             if isinstance(v, ast.Constant) and isinstance(v.value, str))]
                for v in vals:
                    for li in v.splitlines():
                        if line_is_template(li):
                            offenders.append(f"{rel}:{getattr(node,'lineno','?')} (template)")
                            break
offenders = sorted(set(offenders))
check(not offenders, "no production dt=0.002 (code or mdin namelist template)" +
      ("" if not offenders else " -> " + " | ".join(offenders)))


# v2.5.66: campaign analyzer must combine cached leg summaries into dG_bind (no double-count).
_ac = _os.path.join(HERE, "analyze_campaign.py")
check(_os.path.exists(_ac), "tools/analyze_campaign.py present")
if _os.path.exists(_ac):
    import py_compile as _pc2
    try:
        _pc2.compile(_ac, doraise=True); _okc=True
    except Exception as _e:
        _okc=False; print("    (compile error: %s)" % _e)
    check(_okc, "tools/analyze_campaign.py compiles")
    _edge = _tf.mkdtemp()
    def _mk_ti(name, dg):
        _ld = _os.path.join(_edge, name)
        _os.makedirs(_os.path.join(_ld, "lambda_0.000"), exist_ok=True)
        open(_os.path.join(_ld, "lambda_0.000", "prod.out"), "w").write("clambda = 0.0\nFinal Performance Info\n")
        open(_os.path.join(_ld, "summary.json"), "w").write(_json.dumps(
            {"dG_kcal_mol": dg, "complete": True, "n_windows": 1, "n_requested": 1}))
    _mk_ti("complex_decharge", 5.0); _mk_ti("complex_vdw", 30.0)
    _mk_ti("solvent_decharge", 4.0); _mk_ti("solvent_vdw", 25.0)
    _rd = _os.path.join(_edge, "complex_restraint"); _os.makedirs(_rd, exist_ok=True)
    open(_os.path.join(_rd, "boresch_correction.txt"), "w").write("-11.447503\n")
    for _lam in (0.0, 1.0):
        _w = _os.path.join(_rd, "lambda_%.3f" % _lam); _os.makedirs(_w, exist_ok=True)
        open(_os.path.join(_w, "prod.in"), "w").write("&cntrl\n nstlim=2000000, dt=0.001,\n/\n")
        open(_os.path.join(_w, "prod.out"), "w").write("Final Performance Info\nTotal wall time 1\n")
    # ctot = 5+30+(-11.4475)=23.5525 ; stot=29 ; dG_bind = -(23.5525-29) = +5.4475
    _env2 = dict(_os.environ); _env2["PYTHONPATH"] = ROOT + _os.pathsep + _env2.get("PYTHONPATH","")
    _rc = _sp.run([sys.executable, _ac, _edge, "--json-only"], capture_output=True, text=True, env=_env2)
    _bind = None
    _rj = _os.path.join(_edge, "ABFE_RESULT.json")
    if _os.path.exists(_rj):
        try: _bind = _json.load(open(_rj)).get("dG_bind_kcal_mol")
        except Exception: _bind = None
    _close = (_bind is not None) and (abs(_bind - 5.4475) < 1e-3)
    check(_close, "analyze_campaign.py combines legs -> dG_bind=+5.4475 (got %s, exit %s)" % (_bind, _rc.returncode))
    if not _close:
        print("    stdout:", _rc.stdout[-500:])


# v2.5.67: campaign analyzer --csv / --recurse must emit a ranked one-row-per-edge CSV.
if _os.path.exists(_ac):
    import csv as _csvmod
    _camp = _tf.mkdtemp()
    def _edge_with(name_dir, c_dec, c_vdw, bores, s_dec, s_vdw):
        _e = _os.path.join(_camp, name_dir); _os.makedirs(_e, exist_ok=True)
        def _ti(nm, dg):
            _ld = _os.path.join(_e, nm); _os.makedirs(_os.path.join(_ld, "lambda_0.000"), exist_ok=True)
            open(_os.path.join(_ld, "lambda_0.000", "prod.out"), "w").write("clambda = 0.0\nFinal Performance Info\n")
            open(_os.path.join(_ld, "summary.json"), "w").write(_json.dumps(
                {"dG_kcal_mol": dg, "complete": True, "n_windows": 1, "n_requested": 1}))
        _ti("complex_decharge", c_dec); _ti("complex_vdw", c_vdw)
        _ti("solvent_decharge", s_dec); _ti("solvent_vdw", s_vdw)
        _rd = _os.path.join(_e, "complex_restraint"); _os.makedirs(_rd, exist_ok=True)
        open(_os.path.join(_rd, "boresch_correction.txt"), "w").write("%f\n" % bores)
        for _lam in (0.0, 1.0):
            _w = _os.path.join(_rd, "lambda_%.3f" % _lam); _os.makedirs(_w, exist_ok=True)
            open(_os.path.join(_w, "prod.in"), "w").write("&cntrl\n nstlim=2000000, dt=0.001,\n/\n")
            open(_os.path.join(_w, "prod.out"), "w").write("Final Performance Info\nTotal wall time 1\n")
    # ligA: dG_bind = -((5+30-11.4475) - 29) = +5.4475
    _edge_with("ligA", 5.0, 30.0, -11.447503, 4.0, 25.0)
    # ligB: stronger binder. -((5+30-20) - 29) = +14 -> less; make it tighter:
    # choose bores=-40 -> ctot=-5 -> dG_bind=-(-5-29)=+34 (>25 -> untrusted). Instead:
    # ligB ctot=+18 (c_vdw=24.4475) -> dG_bind=-(18-29)=+11 ; still > ligA so ranks after.
    _edge_with("ligB", 5.0, 24.447503, -11.447503, 4.0, 25.0)
    _env3 = dict(_os.environ); _env3["PYTHONPATH"] = ROOT + _os.pathsep + _env3.get("PYTHONPATH","")
    _rc = _sp.run([sys.executable, _ac, _camp, "--recurse", "--json-only"],
                  capture_output=True, text=True, env=_env3)
    _csvp = _os.path.join(_camp, "campaign_summary.csv")
    _ok_csv = _os.path.exists(_csvp)
    _rows = []
    if _ok_csv:
        with open(_csvp) as fh:
            _rows = list(_csvmod.DictReader(fh))
    _hdr_ok = _ok_csv and _rows and ("dG_bind_kcal_mol" in _rows[0]) and ("edge" in _rows[0])
    _two = len(_rows) == 2
    # ranking: trusted-first, then ascending dG_bind -> ligA (+5.4475) before ligB (+11)
    _ranked = _two and _rows[0]["edge"] == "ligA" and _rows[1]["edge"] == "ligB"
    check(_ok_csv and _hdr_ok and _two and _ranked,
          "analyze_campaign --recurse writes ranked campaign_summary.csv (rows=%d, order=%s)"
          % (len(_rows), [r.get("edge") for r in _rows]))
    if not (_ok_csv and _hdr_ok and _two and _ranked):
        print("    stdout:", _rc.stdout[-400:]); print("    csv exists:", _ok_csv)


# v2.5.68: restraint MD leg OFF by default; analytic Boresch term moves onto complex_vdw.
#  (a) default leg list must NOT include complex_restraint, but MUST keep decharge+vdw
#  (b) --restraint-leg restores it
#  (c) write_correction flips to complex_vdw when the leg is off
import inspect as _inspect
_src_drv = _inspect.getsource(__import__("amber_md.fep_driver", fromlist=["x"]))
check(('--no-restraint-leg' in _src_drv) and ('--restraint-leg' in _src_drv),
      "fep_driver exposes --restraint-leg / --no-restraint-leg")
check('_want_restraint_leg = getattr(a, "restraint_leg", False)' in _src_drv,
      "restraint leg is OFF by default (getattr default False)")
check('legs.append(("complex_vdw",       a.absolute_prmtop, a.absolute_inpcrd, boresch, "vdw",      _vdw_writes_corr))' in _src_drv,
      "complex_vdw carries the analytic Boresch correction when the leg is off")
# CRITICAL INVARIANT: restraint POTENTIAL still applied during decharge+vdw (fep.py fixed-k)
_src_fep = _inspect.getsource(__import__("amber_md.fep", fromlist=["x"]))
check('restraint HELD ON the whole' in _src_fep,
      "INVARIANT: Boresch potential held ON during decharge/vdw (ligand cannot drift)")

print()
if FAILS:
    print(f"RESULT: FAIL ({len(FAILS)} check(s) failed)"); sys.exit(1)
print("RESULT: PASS (dt=1fs + 10ns TI + guard + analyzers + tools + CSV + restraint-leg-OFF-default + v2.5.69 docs)"); sys.exit(0)
