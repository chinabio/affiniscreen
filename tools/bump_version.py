#!/usr/bin/env python3
"""bump_version.py NEWVERSION -- update the version across ALL touchpoints at once,
then verify with check_version_sync. Does NOT write the changelog/What's-new prose
(that is intentional human content) but REMINDS you to.

Touchpoints written: amber_md/__init__.py (__version__), VERSION, run_amber.py,
README.md banner. It also checks that Home.py already has a "What's new in vNEW".
"""
from __future__ import annotations
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def sub_file(p: Path, pattern: str, repl: str) -> bool:
    t = p.read_text(); n = re.subn(pattern, repl, t, count=1)
    if n[1]:
        p.write_text(n[0]); return True
    return False

def main(argv) -> int:
    if len(argv) != 2:
        print("usage: bump_version.py NEWVERSION"); return 2
    new = argv[1]
    init = ROOT / "amber_md" / "__init__.py"
    cur = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init.read_text()).group(1)
    print(f"bumping {cur} -> {new}")
    sub_file(init, r'__version__\s*=\s*[\'"][^\'"]+[\'"]', f'__version__ = "{new}"')
    (ROOT / "VERSION").write_text(new + "\n")
    sub_file(ROOT / "run_amber.py", r'__version__\s*=\s*[\'"][^\'"]+[\'"]', f'__version__ = "{new}"')
    sub_file(ROOT / "README.md", r'> \*\*v[0-9][^ ]* \(', f'> **v{new} (')
    hm = (ROOT / "amber_md" / "gui" / "Home.py").read_text()
    if f"What's new in v{new}" not in hm:
        print(f"  REMINDER: add a '#### What's new in v{new}' section to amber_md/gui/Home.py")
        print(f"  REMINDER: add a Changelog entry to amber_md/__init__.py")
        print(f"  REMINDER: add CHANGES_v{new}_*.md")
    print("done. Now run: python tools/check_version_sync.py")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
