# -*- coding: utf-8 -*-
# Part of AffiniScreen.
"""Workflow routing matrix (final62).

Single source of truth for "which page do I use?" across method x engine x
molecule-count. Rendered as a help expander on Results-Compare (and reusable on
any page). Streamlit is imported lazily so the table can be tested headless.
"""
from __future__ import annotations

# (method, engine, count, plan, submit, monitor, results)
ROUTING_ROWS = [
    ("MM-GBSA", "Amber",  "single",
     "Setup & Launch", "Setup & Launch", "Job Monitor", "Results - Single"),
    ("MM-GBSA", "Amber",  "multiple",
     "Setup & Launch", "Setup & Launch  *or*  MM-GBSA Screen",
     "Job Monitor", "MM-GBSA Screen (aggregate) / Results - Compare"),
    ("MM-GBSA", "OpenMM", "single",
     "Setup & Launch", "Setup & Launch", "Job Monitor", "Results - Single"),
    ("MM-GBSA", "OpenMM", "multiple",
     "Setup & Launch", "Setup & Launch  (fan-out)",
     "Job Monitor", "MM-GBSA Screen (aggregate) / Results - Compare"),

    ("ABFE", "OpenFE", "single",
     "Setup & Launch (plan)", "FEP Campaign", "Job Monitor / FEP Campaign",
     "Results - Single / Results - Compare"),
    ("ABFE", "OpenFE", "multiple",
     "Setup & Launch (plan)", "FEP Campaign", "Job Monitor / FEP Campaign",
     "Results - Compare"),

    ("RBFE", "OpenFE", "multiple",
     "Setup & Launch (plan)", "FEP Campaign",
     "Job Monitor / FEP Campaign",
     "FEP Campaign (solve network -> per-ligand dG)"),
]

ROUTING_COLUMNS = ["Method", "Engine", "Molecules",
                   "Plan", "Submit", "Monitor", "Results / Rank"]

NOTES = [
    "RBFE is inherently multi-molecule (a perturbation network); there is no "
    "single-molecule RBFE row.",
    "MM-GBSA submission is engine-asymmetric: OpenMM screens submit ONLY from "
    "Setup & Launch; Amber screens submit from Setup & Launch OR MM-GBSA "
    "Screen. BOTH engines aggregate/rank on MM-GBSA Screen or Results-Compare.",
    "ABFE and RBFE in the GUI are OpenFE-only. FEP Campaign is for OpenFE "
    "NETWORKS (RBFE) and OpenFE ABFE run/gather/solve.",
    "Network RBFE ranking (per-ligand dG via network solve + cycle-closure) "
    "lives on FEP Campaign (OpenFE). Results-Compare ranks PER-LIGAND outputs "
    "(MM-GBSA both engines, OpenFE ABFE) and can also solve+rank a detected "
    "RBFE network inline.",
]


def routing_dataframe():
    import pandas as pd
    return pd.DataFrame(ROUTING_ROWS, columns=ROUTING_COLUMNS)


def render_routing_help(expanded: bool = False):
    import streamlit as st
    with st.expander("Workflow routing matrix - which page do I use?",
                     expanded=expanded):
        st.caption("Method x Engine x molecule-count -> the page for each "
                   "stage. *or* means either page works.")
        st.dataframe(routing_dataframe(), hide_index=True, width="stretch")
        st.markdown("**Notes**")
        for n in NOTES:
            st.markdown(f"- {n}")
