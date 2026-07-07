#!/usr/bin/env python3
"""check_version_sync.py -- fail if the version is not identical across every
touchpoint. Run in CI / pre-release. Source of truth: amber_md/__init__.__version__.

Touchpoints checked:
  * amber_md/__init__.py   __version__
  * VERSION                (whole file, stripped)
  * run_amber.py           __version__
  * README.md              leading banner  > **vX (...)**
  * amber_md/gui/Home.py   a "#### What's new in vX" heading must exist for X
Exit 0 = all in sync; non-zero = drift (prints the offenders).
"""
from __future__ import annotations
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def _truth() -> str:
    t = (ROOT / "amber_md" / "__init__.py").read_text()
    m = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', t)
    if not m:
        print("FATAL: no __version__ in amber_md/__init__.py"); sys.exit(2)
    return m.group(1)

def main() -> int:
    v = _truth()
    bad = []
    # VERSION
    vf = (ROOT / "VERSION").read_text().strip()
    if vf != v: bad.append(f"VERSION = {vf!r} (expected {v!r})")
    # run_amber.py
    ra = (ROOT / "run_amber.py").read_text()
    m = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', ra)
    if not m or m.group(1) != v:
        bad.append(f"run_amber.py __version__ = {(m.group(1) if m else None)!r} (expected {v!r})")
    # README banner
    rd = (ROOT / "README.md").read_text().splitlines()[0]
    if f"v{v}" not in rd:
        bad.append(f"README.md banner does not mention v{v}: {rd!r}")
    # Home.py What's new
    hm = (ROOT / "amber_md" / "gui" / "Home.py").read_text()
    if f"What's new in v{v}" not in hm:
        bad.append(f"Home.py has no \"What's new in v{v}\" section")
    # v2.5.46: stale build artifacts cause `pip install -e .` to report an OLD
    # version (setuptools reuses cached egg-info/PKG-INFO). Fail if any are
    # committed to the source tree, and fail if PKG-INFO (when present) drifts.
    pkg = ROOT / "amber_md.egg-info" / "PKG-INFO"
    if pkg.exists():
        mt = re.search(r'^Version:\s*(.+)$', pkg.read_text(), re.M)
        pv = mt.group(1).strip() if mt else None
        # accept only an exact or pep440-projected match of the truth
        ok = pv in (v, v.replace("a",".post1").replace("b",".post2"))
        if not ok:
            bad.append(f"amber_md.egg-info/PKG-INFO Version = {pv!r} (expected {v!r}); "
                       f"DELETE the egg-info dir -- it shadows the dynamic version")
    stray_pyc = list(ROOT.rglob("*.pyc"))
    if stray_pyc:
        bad.append(f"{len(stray_pyc)} committed .pyc file(s) found (e.g. "
                   f"{stray_pyc[0].relative_to(ROOT)}); remove all __pycache__ dirs")
    if bad:
        print(f"VERSION DRIFT (truth = {v}):")
        for b in bad: print("  -", b)
        return 1
    print(f"version sync OK: every touchpoint = {v}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
