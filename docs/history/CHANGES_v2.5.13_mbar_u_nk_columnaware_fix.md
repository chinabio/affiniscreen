# v2.5.13 - column-aware MBAR u_nk sanitation

## Symptom
MBAR/BAR diverged (dG = -710 vs TI = +8.9 kcal/mol) and pymbar hit DLASCL/SVD
errors even on COMPLETE legs, after the v2.5.12 overflow fix. The headline was
stamped UNTRUSTED (TI/BAR/MBAR disagree).

## Root cause (real data: lig_12944901 solvent_vdw/lambda_0.900/prod.out)
```
MBAR Energy analysis:
  Energy at 0.0000 = ****************   <- Fortran field overflow (decoupled end state)
  Energy at 0.0500 =    -72973.925357   <- clean
  ... every OTHER state clean and smooth ...
```
The soft-core decoupled end-state column is corrupt in ~90% of frames; all
other state columns are fine. The v2.5.12 sanitizer dropped the entire frame
if ANY state cell was bad, so it dropped nearly every frame, emptied ~half the
lambda columns, and left pymbar a degenerate (rectangular/rank-deficient)
reduced-potential matrix. TI was unaffected because it uses dHdl.

## Fix
`FEPAnalyzer._sanitize_u_nk` is now COLUMN-AWARE:
  (a) flag bad cells (non-finite or |u_nk| > _UNK_SANE_MAX_KT kT);
  (b) DROP any STATE COLUMN bad in > _UNK_COL_BAD_FRAC (0.5) of frames --
      unrecoverable, excluded from MBAR/BAR; never median-filled (fabricating
      samples biases dG ~3.7 kcal/mol, also validated);
  (c) drop residual rows with a bad cell in a kept column;
  (d) restrict sampled lambda groups to kept states so the matrix is square.

New constant: `_UNK_COL_BAD_FRAC = 0.5`. Removed: `_UNK_MIN_KEEP` (frame-keep
floor from the superseded drop-frame approach).

## Validation
Synthetic harmonic lambda ladder with a 90%-corrupted end column:
  MBAR(new-sanitized) == MBAR(clean subset) to 0.00000 kcal/mol;
  a clean leg passes through untouched (0 cells flagged, all states kept).

## Operational note
Analysis-only change. Re-run `fep_driver ... --resume --analyze`; no MD re-run
is required for legs whose windows already completed. (Any genuinely crashed
MD windows, e.g. complex_vdw lambda 0.70/0.825, are healed by --resume as usual.)
