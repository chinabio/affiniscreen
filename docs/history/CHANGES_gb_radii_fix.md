# Fix: GB radii not set on parmed-split prmtops (igb=8 / mmpbsa_py_energy fail)

## Symptom
After the numpy.compat fix let parmed load again, the parmed split path ran and
MMPBSA.py then failed:
    CalcError: .../mmpbsa_py_energy failed with prmtop .../topo/complex.prmtop!

## Cause
TopologySplitter._split_parmed built complex/receptor/ligand prmtops with
parmed.strip() but NEVER set GB radii -- so they kept tleap's default 'mbondi'.
mmgbsa.in uses igb=8 (GBn2), which REQUIRES mbondi3 radii. The radii<->igb
mismatch makes mmpbsa_py_energy fail on the complex prmtop. (The ante-MMPBSA.py
path passed --radii=mbondi3 and would have been fine; only the parmed path was
broken.)

## Fix
* _split_parmed now applies parmed.tools.changeRadii(<set>) to every output
  before saving.
* Radius set is chosen from igb (score_mmgbsa):
    igb 1 -> mbondi, 2/5 -> mbondi2, 7 -> bondi, 8 -> mbondi3.
* ante-MMPBSA.py path uses the same `--radii=<set>` (no longer hardcoded).
* score_mmgbsa now calls split(force=True) so a previously-built topo/ with the
  WRONG radii is rebuilt instead of being reused by the idempotent skip.
* MMGBSAAnalyzer.run() now dumps the tails of leftover _MMPBSA_*.out temp files
  on failure, so future failures are self-diagnosing rather than a generic
  CalcError.

## Verified
* igb->radii map matches AmberTools recommendations (mbondi/mbondi2/bondi/mbondi3).
* All mapped sets are valid parmed/ante radius keywords.
* (parmed not installable in the build sandbox; logic verified statically.)
