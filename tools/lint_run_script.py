#!/usr/bin/env python3
"""lint_run_script.py -- build-time linter for generated LSF/bash run scripts.

R1 errexit+pipefail; R2 no unguarded bare command-substitution under errexit;
R3 box-drift escalation markers (warn); R4 non-inverted #BSUB array ranges.
Exit 0 clean, 1 on >=1 ERROR. Usage: lint_run_script.py <script.lsf> [...] | --self-test
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations
import re, sys
ERROR, WARN = "ERROR", "WARN"
def lint_text(text, name="<script>"):
    findings=[]; L=text.splitlines()
    has_errexit=any(s.strip() in ("set -e","set -euo pipefail") or s.strip().startswith("set -eu") for s in L)
    has_pipefail=("set -o pipefail" in text) or ("set -euo pipefail" in text)
    if not has_errexit: findings.append((ERROR,0,"R1: no errexit (set -e / set -euo pipefail)"))
    if not has_pipefail: findings.append((ERROR,0,"R1: no pipefail (set -o pipefail)"))
    errexit=False; bare=re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=\$\(")
    for i,raw in enumerate(L):
        s=raw.strip()
        if s in ("set -euo pipefail","set -e") or s.startswith("set -eu"): errexit=True
        elif s=="set +e": errexit=False
        if not (errexit and bare.match(s)): continue
        if "||" in s or "$((" in s: continue
        if s.startswith(("local ","export ","declare ")): continue
        inner=s.split("$(",1)[1]
        if inner.startswith("set +o pipefail") or inner.startswith("set +e"): continue
        findings.append((ERROR,i+1,f"R2: unguarded bare command-substitution under errexit: {s[:70]}"))
    if "pmemd.cuda" in text and "prod" in text:
        for tag in ("orig_eq.rst","_drift_count","GPU regrid tier"):
            if tag not in text: findings.append((WARN,0,f"R3: escalation marker missing: {tag!r}"))
    for i,raw in enumerate(L):
        m=re.search(r"#BSUB\s+-J\s+\S+\[(\d+)-(\d+)\]",raw)
        if m and int(m.group(1))>int(m.group(2)): findings.append((ERROR,i+1,f"R4: inverted array range [{m.group(1)}-{m.group(2)}]"))
    return findings
def lint_file(path):
    with open(path,encoding="utf-8") as fh: return lint_text(fh.read(),path)
_GOOD="#!/bin/bash\nset -euo pipefail\nmaxTh=$(set +o pipefail; set +e; grep X f | awk '{print}'; true)\n#BSUB -J fep_x[1-16]\npmemd.cuda -O -i prod.in\n# orig_eq.rst _drift_count GPU regrid tier\n"
_BAD="#!/bin/bash\nset -e\nVAR=$(grep X file | head -1)\n#BSUB -J fep_x[16-1]\n"
def _self_test():
    g=[f for f in lint_text(_GOOD) if f[0]==ERROR]; b=[f for f in lint_text(_BAD) if f[0]==ERROR]
    ok=(not g) and len(b)>=2
    print("GOOD errors:",g); print("BAD errors:",b); print("SELF-TEST:","PASS" if ok else "FAIL")
    return 0 if ok else 1
def main(argv):
    if "--self-test" in argv: return _self_test()
    paths=[a for a in argv[1:] if not a.startswith("-")]
    if not paths: print(__doc__); return 2
    rc=0
    for p in paths:
        f=lint_file(p); errs=[x for x in f if x[0]==ERROR]
        for lvl,ln,msg in f: print(f"{lvl} {p}{(':'+str(ln)) if ln else ''}: {msg}")
        if errs: rc=1
        else: print(f"OK {p}: clean ({len(f)} warning(s))")
    return rc
if __name__=="__main__": raise SystemExit(main(sys.argv))
