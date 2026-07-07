# Installation

## 1. Python environment

### Conda (recommended, especially on HPC)
```bash
conda env create -f environment.yml     # creates the 'amber-md' env
conda activate amber-md
pip install -e .[gui,analysis]
```

### pip / venv
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .[gui,analysis]
```

Verify:
```bash
python -c "import amber_md; print(amber_md.__version__)"
```

## 2. Optional extras

| Extra | Install | Needed for |
|-------|---------|-----------|
| `gui` | `pip install -e .[gui]` | Streamlit app (`streamlit`, `nglview`). |
| `analysis` | `pip install -e .[analysis]` | `MDAnalysis`, `pymbar`, `alchemlyb`. |
| `mdtraj` | `pip install "mdtraj<1.10"` | OpenMM MM-GBSA / visualization. Needs numpy≥2, which conflicts with the `numpy<2` Amber pin — install in a separate env if required. |

## 3. Amber & OpenFE

- **Amber / AmberTools** must be available for the MM-GBSA (Amber) path.
  Amber supplies `parmed`/`pytraj` via `PYTHONPATH`; do **not** pip-install them.
- **OpenFE** (and OpenMM + AmberTools in the same env) is required for the
  ABFE/RBFE OpenFE workflows. The GUI runs a **preflight** before launching and
  reports exactly what is missing.

## 4. Cluster (LSF/GPU)

The GUI submits GPU jobs via `bsub`. Set your queue, walltime, and project in
the Setup & Launch **HPC / scheduler** panel. `activate_amber_md.sh` and
`run_gui.sh` help start the app inside the right environment on a login node.

## 5. Version sync

The running version is shown on the Home page and stamped into every generated
`.lsf`. If a job hits an old-format error while the banner shows this version,
the compute/login environment is importing a **stale** `amber_md` — reinstall
this package (`pip install -e .`) into that environment.
