# Changelog

## 2.6.0 — GUI streamlined
- Removed **Amber ABFE** and **Amber RBFE** from the Streamlit GUI:
  - `Setup & Launch`: both combinations disabled in the `COMPAT` matrix; their
    command builders now raise a clear error if reached.
  - `FEP Campaign`: removed the "Amber TI" engine selector — the page is
    OpenFE-only.
  - `Results — Compare`: removed the "Promote MM-GBSA → Amber FEP" action.
  - `Home` / `routing_help`: updated to reflect the supported matrix.
- **Retained** all underlying engine source (`fep_driver`, `rbfe_map`,
  `abfe_*`) for programmatic / CLI use (`amber-fep-driver`).
- Repository reorganized for GitHub: MIT `LICENSE`, `docs/` (user manual,
  install, architecture), `examples/`, and `docs/history/` for the historical
  per-release notes.

## ≤ 2.5.77
See per-release notes in [`docs/history/`](docs/history/).
