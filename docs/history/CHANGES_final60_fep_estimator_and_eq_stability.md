# final60 — two source bugs fixed: MBAR/BAR shape collapse + mid-λ eq blow-up

Applies to `amber_md/fep.py` and `amber_md/config.py` (base: v2.5.4_final60).
Both patched files compile cleanly (`py_compile`) and the Bug 2 fix is
validated against a deliberately-corrupt `u_nk` (one fully non-finite λ window)
through real `alchemlyb` MBAR/BAR.

---

## BUG 2 (pipeline-wide) — `_sanitize_u_nk` corrupted the MBAR matrix

**Symptom in the run:** every leg fell back to TI-only; MBAR/BAR raised
`Shape of passed values is (n,n), indices imply (n+1,n+1)` (and, depending on
the parser path, `SVD did not converge`).

**Root cause** — `fep.py` original lines 1132–1144:

```python
bad = (~np.isfinite(vals)).any(axis=1) | (np.abs(vals-baseline) > _UNK_SANE_MAX_KT).any(axis=1)
return (u_nk if n==0 else u_nk.loc[~bad]), n   # <-- drops whole ROWS (frames)
```

`u_nk` rows are sampled frames, columns are λ-state evaluations. When a
decoupled high-λ window emits soft-core overflow in (almost) every frame
("2499/2500 MBAR energies > 0", non-finite `****`), this filter deletes the
**entire λ group**. After `_decorrelate_per_lambda` concatenates, the matrix has
N−1 populated λ groups but still N state columns (because `mbar_lambda` listed
all N) → rectangular → the alchemlyb shape error. TI survives because its dV/dλ
stream is sanitized separately by `_sanitize_dHdl`.

**Fix (two parts):**

1. `_sanitize_u_nk` now **clamps** offending cells in place (non-finite → row
   baseline; outliers → median ± `_UNK_SANE_MAX_KT`) and **keeps every row**, so
   the matrix stays the same shape and no λ group disappears. Returns
   `(u_nk, n_cells_clamped)`.
2. New `_reconcile_u_nk_states` guard: if a window *genuinely* contributes zero
   surviving frames, restrict **both** index λ-groups and state columns to their
   intersection, so MBAR/BAR always receive a **square (k×k)** matrix. Called
   right before the estimator-fit loop, with a `log.warning` naming the excluded
   states.

**Validation result (this thread):**

| sanitizer | rows dropped | matrix | BAR | MBAR |
|---|---|---|---|---|
| original (drop rows) | 50 (whole λ=0.50) | 4 states × 5 cols | — | **FAILED (SVD did not converge)** |
| patched (clamp+reconcile) | 0 | 5 × 5 (square) | OK | **OK** |

---

## BUG 1 (complex_vdw only) — λ=0.50 eq exploded (T=14,963 K, BOND=3.8M)

**Root cause** — `_eq_in` (`fep.py` original lines 318–335) restarted off the
*gentle* `dens.rst` but reverted to the aggressive integrator:
`dt={dt_ps}` (2 fs), MC barostat (`ntp=1`, default barostat), `tempi=100.0`,
`gamma_ln=2.0`. On a half-decoupled soft-core ligand at mid-λ that regime
diverges. Compounding it, `run_stage` (line 484) only checks the **exit code** —
an eq that integrates to 15,000 K but exits 0 silently passes, and (unlike prod)
eq had no retry/guard.

**Fix (two parts):**

1. `_eq_in` now inherits the density-stage stability via new config knobs:
   `eq_dt_ps=0.001`, `eq_gamma_ln=5.0`, `eq_taup=5.0`, `eq_barostat=1`
   (Berendsen), and restarts at `tempi=temp0` (not 100 K). Step count is scaled
   (`nstlim_eq * dt_ps/eq_dt`) to preserve physical equilibration time.
   Production is unchanged (still MC barostat / `dt_ps`).
2. A **stability gate** after the eq `run_stage` in the LSF array script parses
   the peak `TEMP(K)` from `eq.out` and fails the window (exit 71) if it exceeds
   `EQ_TEMP_MAX_K` (default 1000 K) — so a blown-up box never reaches production.

New `config.py` knobs (safe defaults, backward-compatible via `getattr`):
`eq_dt_ps`, `eq_gamma_ln`, `eq_taup`, `eq_barostat`, `eq_temp_max_K`.

---

## Corrections to the earlier diagnosis (verified against source)

- Line numbers were off by ~1 (`_sanitize_u_nk` starts at **1133**, not 1132).
- `_eq_in` restarts with `irest=1, ntx=5` from `dens.rst` (not from a 100 K
  heat); the destabilising factors are the **2 fs step + MC barostat +
  gamma_ln=2 + tempi=100**, not a cold restart.
- **No "denser vdW schedule" was needed or added.** `vdw_lambdas` already
  includes 0.45/0.50/0.55/0.60/0.65 (config lines 150–152). Adding windows would
  not have fixed λ=0.50 — the integrator did. That proposed change was dropped.

## How to apply & verify
1. Drop in the patched `amber_md/fep.py` and `amber_md/config.py`.
2. Run `run_abfe_smoketest.sh`; confirm MBAR/BAR now report on the healthy legs
   (decharge, solvent_vdw) and that complex_vdw λ=0.50 eq stays < 1000 K.

**Caveat:** clamping restores a *runnable* MBAR, but if a window's cross-energies
are genuinely garbage (real overflow), MBAR will still legitimately disagree with
TI and be rejected by the `_ESTIMATOR_CONSISTENCY_KCAL` gate. Bug 1's fix is what
makes complex_vdw physically sound enough for MBAR to actually agree — you need
both.
