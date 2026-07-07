# v2.5.0 build final38: MM-GBSA complex-leg failure ROOT CAUSE + fix

## Root cause (confirmed on the user's real topology)
MMPBSA failed on the COMPLEX GB calc:
  CalcError: mmpbsa_py_energy failed with prmtop complex.prmtop
Diagnosis from the actual prmtop:
  RADIUS_SET = "modified Bondi radii (mbondi)"   <-- WRONG
  igb in mmgbsa.in = 8 (GBn2), which REQUIRES mbondi3
  IFBOX = 0 (box correctly cleared -- NOT the box bug)
  complex(13394) = receptor(13357) + ligand(37)  (split correct)
So: the prmtop shipped with tleap's default mbondi while igb=8 needs mbondi3.
mmpbsa_py_energy rejects that mismatch -> complex leg aborts on all MPI ranks.

## Why it happened
TopologySplitter._split_parmed applied changeRadii(mbondi3) inside a
try/except that ONLY warned. The radius change silently did not persist
(parmed Structure.radius_set is reliable only on AmberParm), so the prmtop
saved as mbondi and nothing caught it. The ante-MMPBSA fallback path
(--radii) was likewise unverified. (Predates the resname change; unrelated.)

## Fix (amber_md/topology.py)
* _split_parmed: changeRadii failure now RAISES (refuses to save a broken
  prmtop) instead of warning.
* New TopologySplitter._verify_radii(out): reads %FLAG RADIUS_SET directly
  from each saved complex/receptor/ligand prmtop (text-based; works for BOTH
  the parmed and ante-MMPBSA paths) and RAISES on mismatch with the requested
  set (e.g. mbondi vs required mbondi3).
* _verify_radii is called: after the ante-MMPBSA split, and on the cached
  short-circuit in split() (so a previously-broken cache is caught on re-run).

## Verified
  Against prmtop text identical to the user's file:
    Case A (mbondi, igb=8)  -> correctly RAISES with a clear message
    Case B (mbondi3, igb=8) -> PASSES (all three verified)
  All .py compile.

## Immediate workaround for the EXISTING failed run (no re-prep):
  cd ~/Run_dir/run_v250/mmgbsa_20260607_131542
  for f in complex receptor ligand; do
    parmed -p topo/$f.prmtop <<EOF
changeRadii mbondi3
parmout topo/${f}_mbondi3.prmtop
EOF
  done
  # then re-run MMPBSA.py with -cp/-rp/-lp pointing at the *_mbondi3.prmtop files.
