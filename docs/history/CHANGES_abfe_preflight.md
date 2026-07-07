# v2.5.0 build final27: Amber ABFEP environment preflight

## New: amber_md/preflight_abfe.py
A definitive "will Amber ABFEP run here?" check. Run it ON A GPU NODE in the
exact env the FEP job uses (after `module load amber/22` + conda activate):

    bsub -q gpu -gpu "num=1" -W 0:15 -Is bash
    module load amber/22
    conda activate <env>
    python -m amber_md.preflight_abfe --work-dir /tmp/abfe_pf

Does NOT burn real FEP hours. Six checks:
  [1] pmemd.cuda + pmemd.cuda.MPI on PATH         (BLOCK if missing)
  [2] tleap present + builds trivial topology     (BLOCK if missing)
  [3] parmed imports cleanly in THIS env          (BLOCK; the numpy.compat trap
      that broke MM-GBSA when the amber module shadowed the env's parmed)
  [4] pmemd.cuda ACCEPTS GTI/softcore mdin        (BLOCK if the build lacks GTI;
      runs a tiny real 2-residue system with icfe/ifsc/scalpha/scbeta/
      timask/scmask/gti_add_sc/ifmbar -- the exact keyword family fep.py emits)
  [5] pmemd.cuda.MPI multi-window (-ng) smoke     (WARN-only; HREMD path)
  [6] alchemlyb + pymbar importable               (WARN-only; trapezoid-TI
      fallback otherwise)

Exit 0 => GO. Non-zero => >=1 BLOCKER. Writes <work-dir>/abfe_preflight_report.json.

## Why this matters (Amber ABFEP vs MM-GBSA)
The two paths have OPPOSITE env needs:
  * MM-GBSA broke BECAUSE `module load amber/22` shadowed the env's parmed
    (numpy.compat). Fix was to DROP the amber module.
  * ABFEP REQUIRES pmemd.cuda / pmemd.cuda.MPI + GTI, which live in the system
    amber module, NOT the conda env. So it must load amber -- but then the
    topology-build parmed can hit the same numpy.compat trap. Check [3] catches
    exactly that conflict before any GPU hours are spent.

## Verified (stubbed sandbox)
  * Bare env: [1][2][3] BLOCK, [4][5] skip, [6] WARN, exit 1, valid JSON.
  * GO env (fake binaries): [1][2][4][5][6] OK.
  * GTI-rejecting pmemd: [4] BLOCK with 'gti' in blockers (catches non-GTI Amber).
The FEP stack itself (fep.py, fep_driver.py, boresch.py, abfe_topology.py) all
compile; this preflight is the missing runtime gate.
