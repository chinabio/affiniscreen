# v2.5.47 - dependency pins (pandas<3, mdtraj extra)

Date: 2026-06-21

Following the v2.5.46 install on the login node, pip pulled pandas 3.0.3 and warned that mdtraj 1.11 / jaxlib need numpy>=2 (conflicting with the numpy<2 pin Amber parmed requires).

Audit findings (read-only): all mdtraj/MDAnalysis/jax imports are LAZY (inside functions); jax is not imported at all; the FEP analyzer uses only pandas APIs that are 3.0-safe (DataFrame, concat, groupby(level=), reset_index). So nothing was broken. These pins are defensive:
  - pandas<3: avoid an unvalidated major-version jump.
  - mdtraj moved to its own [mdtraj] extra pinned <1.10 (last line that supports numpy<2); removed from [analysis] so `pip install -e .[analysis]` no longer creates a numpy conflict.

No functional code changed; core MD (pmemd+parmed) and FEP analysis (alchemlyb+pandas) paths are identical to v2.5.46.
