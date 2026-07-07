# Examples

Minimal, self-contained inputs to smoke-test the GUI and CLI.

> Replace the placeholder structure files with real inputs for your target.
> They are intentionally tiny so the repo stays lightweight.

## Files
- `protein.pdb` — placeholder receptor (replace with your prepared PDB).
- `ligands.sdf` — placeholder multi-record ligand file (replace with yours).

## Try it (GUI)
1. Launch the GUI: `./run_gui.sh`
2. On **Setup & Launch**, choose **MM-GBSA / Amber**, select
   `examples/protein.pdb` and `examples/ligands.sdf`, then **Validate**.

## Try it (CLI)
```bash
run-amber --help
amber-fep-driver --help
```
