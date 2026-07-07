"""AffiniScreen -- Streamlit GUI (multipage)."""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of the AffiniScreen.

from __future__ import annotations
import streamlit as st

# Single source of truth for the product name -- change here to rebrand.
APP_NAME = "AffiniScreen"
APP_TAGLINE = "Binding free-energy screening: MM-GBSA triage \u00b7 FEP (ABFE/RBFE) confirmation"

st.set_page_config(page_title="Home",
                   layout="wide", page_icon="bio",
                   initial_sidebar_state="expanded")

import amber_md

# Sidebar: give the landing page a friendly "Home" label at the top of the nav.
with st.sidebar:
    st.markdown("## \U0001F3E0 Home")
    st.caption(APP_NAME)
    # Always-visible version badge so it is obvious which code is running.
    _build = getattr(amber_md, "__build__", "")
    st.caption(f"**v{amber_md.__version__}**" + (f" ({_build})" if _build else ""))

st.title(APP_NAME)
_build = getattr(amber_md, "__build__", "")
st.caption(f"v{amber_md.__version__}" + (f" ({_build})" if _build else "")
           + f" - {APP_TAGLINE}")

st.markdown("""
### New here? Start in 3 steps

1. **Settings** - point the app at your cluster and tools (do this once).
2. **Setup & Launch** - choose a method, pick your protein + ligands, and launch.
3. **Results - Compare & Rank** - when jobs finish, rank your ligands.

> Not sure which method to pick? See **"Which method should I use?"** below.

---

### Welcome

Pick a workflow from the **sidebar** (listed in run order):

| Page | Use when... |
|------|-------------|
| **Settings** | Configure your environment once -- scheduler (LSF; SLURM planned), GPU/CPU queues, Amber & OpenFE locations, and default paths. Other pages read these as their defaults. Start here when adapting the workflow to a new machine or cluster. |
| **Setup & Launch** | Configure and launch any supported run -- MM-GBSA (Amber or OpenMM), ABFE (OpenFE), or RBFE (OpenFE); single or batch. The unified wizard. |
| **Job Monitor** | Track running / pending scheduler jobs and tail their logs. |
| **MM-GBSA Screen** | Submit / resume / aggregate a multi-ligand MM-GBSA screen (Amber or OpenMM). Use **Results - Compare** to rank it. |
| **FEP Campaign** | Run / gather / solve an OpenFE RBFE/ABFE free-energy network with cycle-closure QC, and **inspect atom mappings**. |
| **Results - Single Molecule** | Detailed energetics and convergence for one finished run. |
| **Results - Compare & Rank** | Rank a campaign of ligands and score vs experiment. |

---

#### Supported workflows

| Method | Engine | Status |
|--------|--------|--------|
| **MM-GBSA** | Amber (AmberTools MMPBSA.py) | Supported |
| **MM-GBSA** | OpenMM + AmberTools | Supported (experimental) |
| **ABFE** | OpenFE (AbsoluteBindingProtocol) | Supported |
| **RBFE** | OpenFE (relative binding network) | Supported |

> Amber ABFE and Amber RBFE are **not** exposed in the GUI. The underlying
> engine code (`amber_md.fep_driver`, `amber_md.rbfe_map`, and the `abfe_*`
> modules) is retained in the package for programmatic / CLI use.

---

#### Which method should I use?

| If you want to... | Use | Speed | Accuracy |
|---|---|---|---|
| Quickly triage many ligands against one target | **MM-GBSA** | Fast | Approximate (ranking) |
| Get an absolute binding free energy for each ligand | **ABFE** | Slow | High |
| Precisely rank a series of similar ligands | **RBFE** | Slow | High (relative) |

- **MM-GBSA** is the cheap first pass - run it on your whole library, then
  promote the top hits to a rigorous free-energy method.
- **RBFE** shines when your ligands are chemically similar (a congeneric
  series); it computes the *difference* between pairs, which cancels errors.
- **ABFE** gives a standalone number per ligand - use it for diverse ligands or
  when you need an absolute value, not just a ranking.
""")

from amber_md.gui.common import get_lsf_jobs
st.divider()
st.subheader("Your jobs")
jobs = get_lsf_jobs()
if jobs:
    import pandas as pd
    df = pd.DataFrame(jobs)
    st.dataframe(df, hide_index=True, width='stretch')
    st.caption(f"{len(jobs)} job(s) active. Refresh page to update.")
else:
    st.info("No active jobs right now. Launch one from **Setup & Launch**.")