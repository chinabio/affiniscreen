#!/usr/bin/env python3
"""validate_mdin.py (amber_md v2.5.31) -- CLI wrapper around amber_md.mdin_validator."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from amber_md.mdin_validator import validate_paths, format_issues  # noqa: E402


def main(argv) -> int:
    roots = argv[1:] or ["."]
    errs, warns, n = validate_paths(roots)
    print(f"validate_mdin: inspected {n} window(s); "
          f"{len(errs)} error(s), {len(warns)} warning(s)\n")
    body = format_issues(errs, warns)
    if body:
        print(body)
    if not errs:
        print("\nPASS: no parser-fatal mdin issues found. Safe to submit." if n
              else "\nNo windows found to validate (check the path).")
    else:
        print(f"\nFAIL: {len(errs)} error(s) would abort pmemd at setup. "
              "Fix before submitting.")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
