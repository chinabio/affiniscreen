# Fix: system amber module shadows env AmberTools on GPU queue

## Symptom
OpenMM MM-GBSA submitted to the LSF gpu queue ran MD fine, then failed at the
topology-split / MMPBSA step:
    ModuleNotFoundError: No module named 'numpy.compat'
from parmed (and the ante-MMPBSA.py fallback hit the same).

## Cause
The LSF #BSUB script did `module load amber/22.8`, which prepends
/share/apps/amber/rhel8/amber22/bin to PATH and makes `parmed` / `ante-MMPBSA.py`
resolve to Amber's OLD system install. That parmed does
`from numpy.compat import asbytes`, removed in NumPy >= 1.24, so it crashes when
a modern numpy is visible. The the login node (--submit local) run worked because the
amber module was NOT loaded there; the env's own AmberTools were used.

## Fix (mmgbsa_openmm._submit_to_lsf)
* Any `amber*` module is DROPPED from the GPU job's module list (CUDA/gcc kept).
  The OpenMM MM-GBSA pipeline is fully self-contained in the conda env (its own
  tleap/antechamber/MMPBSA.py/parmed), which is internally consistent.
* New `--modules` CLI flag overrides the module list entirely; amber* is still
  filtered out.
* Defensive PATH scrub removes any /share/apps/amber*/ or /amber*/bin entries
  that a site profile might re-add after conda activate.

## Verified
* Generated #BSUB script contains gcc + cuda, NO `module load amber`, the PATH
  scrub line, and the inner command forced to `--submit local`.
