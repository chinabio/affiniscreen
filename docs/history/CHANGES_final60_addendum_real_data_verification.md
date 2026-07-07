# final60 — addendum: real-data verification & HREMD gate gap closed

This addendum supplements `CHANGES_final60_fep_estimator_and_eq_stability.md`
after validating the fix against two **real** pmemd `eq.out` files from the
failing run (`complex_vdw` λ=0.45 and λ=0.50) and auditing every code path that
generates an eq stage or consumes `u_nk`.

---

## 1. Real-data reproduction of the exact MBAR shape error

Using the real λ=0.50 `eq.out` (the window that detonated to 14,963 K), parsed
with the genuine `alchemlyb.parsing.amber.extract_u_nk` and run through the
real `decorrelate → sanitize → fit` pipeline:

| | OLD (`u_nk.loc[~bad]`, no reconcile) | PATCHED (clamp + `_reconcile_u_nk_states`) |
|---|---|---|
| corrupt cells | dropped 235 rows | clamped 5,207 cells, **0 rows dropped** |
| matrix | values (3,3) vs index (4,4) — **rectangular** | (442 × 4) — **square** |
| BAR | **FAILED:** `Shape of passed values is (3,3), indices imply (4,4)` | **OK** (dG = −0.108) |
| MBAR | **FAILED:** same error | **OK** (dG = +0.494) |

This reproduces the **exact** `Shape of passed values is (n,n), indices imply
(n+1,n+1)` error from the run, from real data, and shows the patched pipeline
recovers both estimators. Codified as a regression test —
`tests/test_fep_estimator_recovery.py`.

## 2. Corrected understanding of the drop-rows failure

The earlier writeup implied a corrupt window contributes "one bad frame."
Real-data measurement showed worse: a single `****` (→ non-finite) detonation
frame **inflates the per-frame median**, dragging ~60 % of otherwise-healthy
frames outside `_UNK_SANE_MAX_KT`. The old sanitizer therefore dropped
**146–235 rows**, not 1 — collapsing the matrix far more readily than described.
The clamp approach is immune to this because it never lets a non-finite cell
define the baseline (it falls back to the global finite median).

## 3. Parser behaviour confirmed (not assumed)

`extract_u_nk` does **not** raise on a `****` MBAR block — it returns a frame
with non-finite cells (logging "N MBAR energies are > 0.0" and "prematurely
terminated run"). So the corrupt data reaches `_sanitize_u_nk` as NaN/inf, which
is exactly what the clamp path is built to neutralise. `extract_dHdl` returned
all-finite dV/dl here, consistent with TI surviving while MBAR/BAR died.

## 4. NEW FIX — HREMD per-window eq had no stability gate

The original Bug 1 patch added the eq temperature gate only to the
**single-window LSF-array** path (`build_lsf_array`). The **HREMD** path
(`build_lsf_hremd`) runs its own per-window `min → dens → eq` loop *before*
replica exchange, and that loop was **ungated** — a λ=0.50-style blow-up would
exit 0, write a poisoned `eq.rst`, and corrupt every replica in the exchange.

**Fix:** the same peak-`TEMP(K)` gate (`EQ_TEMP_MAX_K`, default 1000 K, exit 71)
is now emitted inside the HREMD per-window eq loop. Verified: **two** eq
stability gates are now present in `fep.py` (LSF-array + HREMD).

## 5. Resume path audited — safe, no change needed

`build_lsf_array_resume` regenerates the script from `build_lsf_array` (which now
contains the gate) and only rewrites the `#BSUB -J ...[range]` line via a
targeted `re.sub`. The gate text is preserved verbatim under `--resume`.

---

## Net file changes in this addendum
- `amber_md/fep.py` — **+ HREMD eq stability gate** (on top of the original
  final60 patches). Compiles clean; 1666 lines.
- `tests/test_fep_estimator_recovery.py` — new regression test (real fixtures).
- `tests/fixtures/eq_complex_vdw_lambda_0.450.out` — healthy eq (peak 300.7 K).
- `tests/fixtures/eq_complex_vdw_lambda_0.500.out` — detonated eq (14,963 K,
  all-`****` MBAR block).

## Run the regression test
```bash
pip install alchemlyb pymbar pandas numpy scipy pytest
pytest -q tests/test_fep_estimator_recovery.py
```
