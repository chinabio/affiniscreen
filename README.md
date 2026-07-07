# AffiniScreen

A Streamlit application and Python package for structure-based binding
free-energy screening. It provides one unified GUI to configure, launch,
monitor, and analyze three classes of calculation on an LSF/GPU cluster:

| Method | Engine | Description |
|--------|--------|-------------|
| **MM-GBSA** | Amber (AmberTools `MMPBSA.py`) | Endpoint binding-energy screening. |
| **MM-GBSA** | OpenMM + AmberTools | OpenMM MD with AmberTools scoring (experimental). |
| **ABFE** | OpenFE (`AbsoluteBindingProtocol`) | Absolute binding free energy. |
| **RBFE** | OpenFE (relative binding network) | Relative binding free-energy network with cycle-closure QC. |

> **Scope note.** Amber ABFE and Amber RBFE are **not** exposed in the GUI.
> The underlying engine code (`amber_md.fep_driver`, `amber_md.rbfe_map`, and
> the `abfe_*` modules) is retained in the package for programmatic / CLI use,
> but the graphical workflows are OpenFE-based for ABFE/RBFE and Amber/OpenMM
> for MM-GBSA. See [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md).

> ⚠️ **Before you run: configure your environment.** AffiniScreen ships
> with **placeholder** scheduler settings. Open the **Settings** page in
> the GUI (or edit `site_config.yaml`) and set your **scheduler project /
> account** (the default is `your-project`), **GPU/CPU queues**, and the
> paths to **Amber** and **OpenFE**. Cluster jobs will be rejected until a
> valid project/queue is set for your site.

---

## Quick start

```bash
# 1. Clone
git clone <your-repo-url> affiniscreen
cd affiniscreen

# 2. Create the environment (conda recommended on HPC)
conda env create -f environment.yml          # creates the 'amber-md' env
conda activate amber-md
pip install -e .[gui,analysis]               # install this package + extras

# 3. Launch the GUI
./run_gui.sh                                  # opens http://localhost:8501
```

Or with pip only (no conda):

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .[gui,analysis]
streamlit run amber_md/gui/Home.py
```

See [`docs/INSTALL.md`](docs/INSTALL.md) for HPC / Amber / OpenFE details.

---

## The GUI at a glance

The app is multipage; the sidebar lists pages in run order:

1. **Setup & Launch** — pick *Method × Engine × Scope*, choose inputs
   (protein + ligands), review a method-specific parameter panel, and launch.
2. **Job Monitor** — track LSF jobs and tail logs.
3. **MM-GBSA Screen** — submit / resume / aggregate a multi-ligand MM-GBSA run.
4. **FEP Campaign** — OpenFE run → gather → solve → cycle-closure, plus
   atom-mapping inspection.
5. **Results — Single Molecule** — detailed energetics and convergence.
6. **Results — Compare & Rank** — rank a campaign and score against experiment.

A full, screenshot-illustrated walkthrough is in
[`docs/USER_MANUAL.md`](docs/USER_MANUAL.md).

---

## Command-line entry points

Installing the package exposes:

```bash
run-amber --help            # run_amber.py: MM-GBSA / MD driver
amber-fep-driver --help     # amber_md.fep_driver: FEP engine (advanced/CLI)
```

The `amber-fep-driver` still supports Amber ABFE/RBFE for advanced users on the
command line; those paths are simply not surfaced in the GUI.

---

## Repository layout

```
affiniscreen/
├── amber_md/               # the Python package
│   ├── gui/                # Streamlit multipage app (Home.py + pages/)
│   ├── fep_driver.py       # FEP engine driver (CLI-supported)
│   ├── rbfe_map.py         # RBFE network planner
│   ├── abfe_*.py           # ABFE integration/restraint/self-heal modules
│   └── ...                 # config, MD, MM-GBSA, batch, analysis helpers
├── docs/                   # documentation (see below)
│   ├── USER_MANUAL.md      # illustrated GUI walkthrough
│   ├── INSTALL.md          # install / environment / HPC setup
│   ├── images/             # screenshots (placeholders to fill in)
│   └── history/            # historical per-release change notes (provenance)
├── examples/               # tiny sample inputs + a runnable example
├── tools/                  # standalone validators / dev utilities
├── tests/                  # test suite
├── environment.yml         # conda environment
├── pyproject.toml          # packaging (pip install -e .)
├── run_gui.sh              # convenience launcher
├── CHANGELOG.md            # release history (summary)
├── CONTRIBUTING.md
└── LICENSE                 # MIT
```

---

## Documentation

| Doc | What it covers |
|-----|----------------|
| [`docs/USER_MANUAL.md`](docs/USER_MANUAL.md) | Step-by-step GUI usage with screenshots. |
| [`docs/INSTALL.md`](docs/INSTALL.md) | Environment, Amber, OpenFE, and cluster setup. |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | How the GUI maps selections to engine commands. |
| [`README_CLEAN_RUN.md`](README_CLEAN_RUN.md) | Pre-submission checklist for FEP legs. |
| [`docs/history/`](docs/history/) | Historical release-by-release change notes. |

---

## License

Released under the [MIT License](LICENSE).

*Version 2.6.0.*
