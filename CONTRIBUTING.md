# Contributing

Thanks for improving the AffiniScreen!

## Development setup
```bash
pip install -e .[gui,analysis]
```

## Ground rules
- **Keep the GUI a thin layer.** Science/engine behavior should live in the
  engine modules, not in Streamlit pages.
- **Don't re-enable Amber ABFE/RBFE in the GUI** without validating them
  end-to-end first (that is why they are GUI-disabled).
- Run the tests before opening a PR:
  ```bash
  python -m pytest -q
  ```
- Compile-check the GUI after edits:
  ```bash
  python -m py_compile amber_md/gui/**/*.py
  ```

## Versioning
`amber_md/__init__.py` is the single source of truth (`__version__`). Keep
`VERSION` and the README in sync (`tools/check_version_sync.py`).
