# Fix: GB calc aborts on periodic complex.prmtop (IFBOX != 0)

## Symptom (revealed by the new _MMPBSA_*.out diagnostic dump)
mmpbsa_py_energy on complex.prmtop failed with:
    Error: gb>0 is incompatible with periodic boundary conditions.
    Error: To use this method set IFBOX in the PRMTOP file to 0.

## Cause
The dry complex/receptor/ligand prmtops are produced by parmed.strip() from the
*solvated* prmtop, which carries a periodic box (IFBOX=1). GB (gb=8) requires a
NON-periodic topology, so the energy program aborts before computing anything.
(The radii fix from the previous build was correct and is confirmed in the log:
"parmed split applied GB radii 'mbondi3'".)

## Fix (topology._split_parmed)
Before saving each dry GB prmtop, set `struct.box = None` so IFBOX=0, in
addition to applying the GB radius set. The solvated.prmtop (-sp reference) is
left periodic, which is correct.

## Verified
* topology.py compiles; box-removal + radii applied to complex/receptor/ligand.
* (parmed not installable in build sandbox; the failing Amber error message
  explicitly prescribes this exact fix -- set IFBOX=0.)
