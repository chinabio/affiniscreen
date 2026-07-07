# amber_md v2.5.70 — complex_decharge protocol fix + MD-only NaN guard

**Base:** v2.5.69

## TL;DR
The `complex_decharge` free-energy analysis was producing garbage
(TI = -56.6 vs BAR = -19.3 kcal/mol; MBAR refused to solve with
`LinAlgError: SVD did not converge`). Root-caused to the **decharge alchemical
protocol**, not the analysis code. Fixed by decharging with **hard van der Waals**
and a **densified high-lambda schedule**, plus a long-standing **false-positive in
the production NaN guard** that discarded healthy high-lambda windows.

The underlying MD was healthy throughout — all 15 windows ran the full
5,000,000 steps with 100% finite trajectories. The failures were thermodynamic
(protocol) + reporting (guard), not simulation crashes.

---

## ROOT CAUSE (complex_decharge)
The decharge production namelist contained:

```
icfe=1, ifsc=1, clambda=0.825,
timask1=':LIG', timask2='',
scmask1=':LIG', scmask2='',          # ligand IN the vdW soft-core region
                                     # (no crgmask on decharge stage)
gti_ele_sc=1, gti_vdw_sc=1,          # vdW softened
gti_add_sc=1, gti_chg_keep=1,        # soft-core charges KEPT...
```

So the ligand was simultaneously (a) having its charges removed by `clambda`
**and** (b) sitting in a **softened** vdW core (`gti_vdw_sc=1`) with
`gti_chg_keep=1`. With vdW softened, atoms may partially overlap; the partial
charges being scaled by `clambda` then approach through that softened core, giving
a `-q*q/r` **charging singularity**. Diagnostic evidence on the real run:

* MBAR table `Energy at 0.000` re-evaluations exploded to `4.8e8` / `****`
  (Fortran field overflow) for windows lambda >= 0.725 (116 -> 5352 overflow cells).
* `<dV/dl>` swung **+102 -> -507 -> -88 kcal/mol** across 0.625 -> 0.95.
* Adjacent-pair BAR work reached **600-700 kT** with `err = 0.000 / nan` and
  sign flips — the signature of zero phase-space overlap.
* alchemlyb MBAR: `LinAlgError: SVD did not converge` (singular u_nk) — hence the
  **missing MBAR row** in `dG_estimators.csv`.

(The same family of overlap failure historically hit `complex_vdw` on the coarse
schedule; that leg already moved to the dense 44-window grid and is being re-run.
v2.5.70 does **not** change the vdw protocol.)

## FIX

### 1. Stage-aware decharge soft-core (`config.py`, `fep.py`)
New config fields, defaulting to the corrected behaviour:

| field | default | effect on DECHARGE stage |
|-------|---------|--------------------------|
| `decharge_gti_vdw_sc`   | `0` | **hard vdW** — ligand vdW NOT softened during decharge |
| `decharge_gti_chg_keep` | `0` | charges genuinely removed (not kept on soft-core) |

`FEPSetup._stage_softcore_block()` is now stage-aware: it emits
`gti_vdw_sc=0, gti_chg_keep=0` for the **decharge** stage and the unchanged
`gti_vdw_sc=1, gti_chg_keep=1` for the **vdw** stage. `ifsc=1` is still forced on
the decharge stage (mismatched single-topology mask atom counts require it), and
electrostatic soft-core (`gti_ele_sc=1`) is retained — only the *vdW* softening
and the contradictory charge-keep are turned off where they caused the singularity.

### 2. Densified decharge schedule (`config.py`)
`decharge_lambdas`: **15 -> 21 windows**, finer across the high-lambda charging
region (0.5 -> 0.975). The legacy grid is preserved as `decharge_lambdas_legacy`
for reproducibility / `--decharge-lambdas` overrides.

```
old: 0.0 0.025 0.05 0.1 0.175 0.275 0.375 0.5 0.625 0.725 0.825 0.9 0.95 0.975 1.0
new: 0.0 0.025 0.05 0.1 0.175 0.275 0.375 0.5 0.5625 0.625 0.6875 0.725
     0.775 0.825 0.8625 0.9 0.925 0.95 0.9625 0.975 1.0
```

### 3. MD-trajectory-only NaN/`****` production guard (`fep.py`, ~line 1207)
The v2.5.59 guard grepped the **entire** `prod.out` for `NaN` / `**********`. But
the `MBAR Energy analysis:` tables legitimately contain `****` overflows and huge
finite cross-terms at distant off-diagonal lambda pairs (handled downstream by
`_sanitize_u_nk`). Those are **not** detonations — yet the guard flagged healthy
high-lambda windows (0.725-1.000) as `rc=0-but-NaN` and discarded good production.
The guard now strips everything from the first `MBAR Energy analysis:` line onward
(`sed '/MBAR Energy analysis:/,$d'`) and checks only the real per-step MD energy
records. A genuine single-step detonation (TEMP=NaN / Etot=NaN / BOND `****` in the
MD section) is still caught and still returns 255.

---

## MIGRATION / ACTION REQUIRED
* **Re-run `complex_decharge`** from scratch with v2.5.70 — both the namelist and
  the lambda set changed, so old windows are not MBAR-compatible.
* `complex_vdw` is unchanged by this release; continue the in-flight 44-window run.
* No API changes. New config fields have safe defaults; existing overrides honoured.

## VERIFICATION (this release)
* Generated decharge `prod.in` now emits `gti_vdw_sc=0, gti_chg_keep=0`,
  21-window `mbar_lambda`, `ifsc=1`, `gti_ele_sc=1` (checked in build).
* Generated vdw `prod.in` is byte-identical to v2.5.69 except version-string
  context (no protocol drift).
* NaN guard: a synthetic prod.out with `****` ONLY in the MBAR table passes;
  a synthetic prod.out with `Etot = NaN` in the MD section still returns 255.
* `tools/check_version_sync.py` passes (all touchpoints = 2.5.70).

## Version track
`amber_md/__init__.py`, `run_amber.py`, `VERSION`, README banners -> **2.5.70**.
