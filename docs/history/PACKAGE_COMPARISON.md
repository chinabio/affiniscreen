# ABFE Reference-Package Comparison

**Document version:** 1.1 (paired with amber_md_workflow v2.5.57)
**Date:** 2026-06-21
**Purpose:** A permanent, grounded reference so the external ABFE codebases do
**not** need to be re-uploaded. Every parameter below was extracted directly
from the source files of each package as uploaded to this project. Lines marked
**[NOT GROUNDED]** could not be verified from an uploaded artifact and must not
be quoted as fact.

This complements `COMPARATIVE_ANALYSIS.md` (which focuses on the FEP+ `.msj`
protocol and the pmemd `skinnb` root-cause). This file is the broad
side-by-side of **all** reference packages.

---

## 0. Packages covered

| Short name | Engine | Source artifact in project | Grounded? |
|---|---|---|---|
| **This workflow** | Amber pmemd.cuda | `amber_md_workflow_v2.5.xx` | ✅ full |
| **FEP-SPell-ABFE** | Amber pmemd (GTI) | `FEP-SPell-ABFE-main.zip` | ✅ full (code read) |
| **the reference platform** | OpenMM | `reference-platform.zip` | ✅ full (code read) |
| **FEP+** | Schrödinger Desmond | `ATO_SEPTTR_*.msj` + `multisim.log` | ✅ full (msj read) |
| **BAT.py** | Amber / OpenMM | `BAT.py-master.zip` (**96-byte stub**) | ❌ NOT uploaded — see §6 |

| **BAT.py** | Amber / OpenMM | `BAT.py-master.zip` (7 MB, full source) | ✅ full (code read) |

## 1. Master comparison table (all grounded unless noted)

| Dimension | This workflow | FEP-SPell-ABFE | the reference platform | FEP+ |
|---|---|---|---|---|
| MD engine | pmemd.cuda | pmemd (GTI/CUDA) | OpenMM | Desmond |
| Alchemy | dual-topology TI + MBAR | dual-topology TI, `ifsc=1`, MBAR | OpenMM softcore, MBAR | `desmond:fep` REST |
| Soft-core | scalpha/scbeta (stage-dep) | **scalpha=0.2, scbeta=50.0** | OpenMM nonbonded softcore | built-in |
| Timestep | 2 fs (`dt=0.002`) | **2 fs** (4 fs if HMR) | **2 fs** (`dt_ps=0.002`) | RESPA `[4,4,8]` fs |
| Temperature | 298 K | **298.15 K** | **298.15 K** (`targetT_K`) | 298 K (10→298 anneal) |
| Cutoff | 10 Å | **9 Å** (`cutoff: 9.`) | OpenMM default | 7–9 Å region |
| Barostat | Berendsen (taup=5, v2.5.54) | **Berendsen taup=2.0** | **Monte-Carlo** (`npt_mc_freq=25`) | MTK → NVT/μVT |
| Production / window | 5 ns | **5 ns** (`complex_length: 5.`) | snapshot-driven (see §4) | **10 ns** (`time=10000`) |
| Replica exchange | optional HREMD | optional (`remd`, 1000 exch) | **REPEX** (`repex.py`) | **REST + λ-hopping** (HREX) |
| Water sampling | standard NPT | standard NPT | standard NPT | **GCMC μVT** |

### λ-window counts (grounded)

| Leg | This workflow | FEP-SPell | the reference platform (recipe) | FEP+ |
|---|---|---|---|---|
| Complex vdW (decouple) | per config | **44** | up to **44** (`v44`) | part of 108 |
| Solvent vdW (decouple) | per config | **44** | recipe-dependent | part of 60 |
| Complex restraint | **21** | **16** | recipe `r01/r02` = **8** | part of 108 |
| Complex electrostatics | per config | (charge corr.) | up to **29** (`e29`) | charge pathway |
| **Total complex** | — | ~88–104 | recipe-dependent | **108** |
| **Total solvent** | — | ~44 | recipe-dependent | **60** |

> Your **21-window restraint schedule** is clearly an expansion of FEP-SPell's
> **16-window** `abfe_complex_restraint_lambdas` — same 0.0→1.0 shape, denser
> near the endpoints. This confirms FEP-SPell as the lineage of this workflow.

---

## 2. FEP-SPell-ABFE (grounded from code)

Files read: `abfe/md/amber_alchemy_mdin.py`, `testing/abfe_testing/config.yaml`,
`abfe/md/amber_mdin.py`.

- **MD control** (`amber_alchemy_mdin.py`): `taup=2.0`, `gamma_ln=2.0`,
  `ifsc=1`, `scalpha=0.2`, `scbeta=50.0`, `ifmbar=1` with `mbar_states` /
  `mbar_lambda`. Restraint-removal stage uses Boresch via `timask1/2` +
  `clambda`. Eq stage `nstlim=10000`.
- **config.yaml**: `temperature=298.15`, `cutoff=9.`, `gti_add_sc=1`
  ("double-decoupled"; `=3` would be double-annihilation), eq+prod legs
  `solvated_length=5` ns and `complex_length=5` ns, `timestep=0.002`
  (0.004 with HMR), optional `remd` with `remd_numexchg=1000`.
- **Boresch force constants**: `bond_const=10.0` kcal/mol/Å², `angle_const=100.0`
  kcal/mol/rad², `dihedral_const=100.0` kcal/mol/rad²; positional-restraint
  `restraint_wt=10` (eq) / `5` (prod) kcal/mol/Å².
- **λ schedules** (config.yaml, counted): solvated-vdw **44**, complex-vdw
  **44**, complex-restraint **16**.
- **Analysis**: `methods: ["MBAR", "TI"]`; explicit `solvent_mask` for
  charge-correction.

**Relationship to this workflow:** same engine family, same dual-topology TI +
MBAR design, same 5 ns legs. This workflow **diverged on the barostat**
(FEP-SPell `taup=2.0` → this workflow `taup=5.0` + `skinnb=3.0`) specifically to
cure the pmemd.cuda fixed-grid halt that FEP-SPell's defaults are also
theoretically exposed to but did not trigger on its smaller test systems
(benzene/phenol in 181L).

---

## 3. the reference platform (grounded from code)

Files read: `reference_platform/protocols/abfe/recipes.py`,
`reference_platform/configs/global_keys.py`, `reference_platform/protocols/abfe/config_types.py`,
`examples/abfe/abfecfg.yaml`, `examples/abfe/run.sh`.

- **Engine**: OpenMM. Repex driver `reference_platform/app/dyn/repex.py`; vanilla MD
  `vanilla.py`. Multi-GPU (`CUDA_VISIBLE_DEVICES=0..7` in `run.sh`).
- **Integrator defaults** (`global_keys.py`): `dt_ps=0.002`,
  `nstep_per_snapshot=500`, `targetT_K=298.15`, **MC barostat**
  `npt_mc_freq=25`.
- **λ is a 3-vector** `[vlam, elam, reslam]` (vdW, electrostatics, restraint) —
  the reference platform composes the full ABFE path from three recipe families in
  `recipes.py`:
  - electrostatics: `e05`(5), `e11`(11), `e23`(23), `e29`(29) states
  - vdW: `v05`(5), `v18`(18), `v20`(20), `v24`(24), `v33`(33), **`v44`(44)** states
  - restraint: `r01`(8), `r02`(8) states
- **Boresch** (`global_keys.py` GKBoresch): default `k_kcal=25.0`,
  `r_theta_phi` + `alpha_beta_gamma` equilibrium values, `k_r_a_dih_kcal` force
  constants — same Boresch 6-DOF scheme this workflow uses.
- **Sampling sizing** (`config_types.py`): `md_nsnapshots=2000`,
  `md_eq_nsnapshots=600` (production length is snapshot-count × 500 steps ×
  2 fs = **2 ns of sampling** at defaults; **[NOT GROUNDED]**: any per-system
  override of these defaults).
- **Stages** (`protocols/abfe/`): `stage_makebox`, `stage_boresch`,
  `stage_sysA`, `stage_sysB`, `stage_post_analysis`; analysis via
  `main_fe_mbar.py` (MBAR) and `main_fe_restraints.py`.

**Why the reference platform never hits the "box changed too much" halt:** OpenMM's MC barostat +
automatic nonbonded regridding means there is no fixed pairlist grid to violate.
This is the structural reason the v2.5.54 `skinnb` fix is Amber-specific.

---

## 4. FEP+ (grounded from `.msj` — summary; see COMPARATIVE_ANALYSIS.md for detail)

- Complex leg **108 λ windows**, solvent leg **60**; production **10 ns/window**
  (`time=10000.0`), RESPA `[4,4,8]` fs.
- **REST + λ-hopping** (HREX, `exchange_probability=0.3`) and **GCMC** water
  (`solvent_density=0.03248`, `mu_excess=-6.137`).
- `scale_solvent_vdw=0.75`. Pre-FEP relaxation ladder 100/12/20/20 ps then 1 ns.
- Full multisim wall time **12h 23′ 56″** for one ligand (`multisim.log`).

---

## 5. Design takeaways

1. **Lineage**: this workflow is an Amber/pmemd implementation in the
   **FEP-SPell-ABFE** tradition (dual-topology TI, `ifsc=1`, scalpha/scbeta,
   MBAR, Boresch, 5 ns legs), with window schedules expanded (restraint 16→21).
2. **Sampling philosophy differs sharply by engine**: FEP+ spends ~10 ns ×108
   windows with REST+GCMC; the Amber/OpenMM open tools use fewer ns but more or
   comparable windows and rely on MBAR/repex. Cost is not directly comparable
   without normalizing windows × ns × GPU.
3. **The box-halt is unique to fixed-grid Amber**: neither the reference platform (OpenMM MC
   barostat) nor FEP+ (Desmond re-tiling) can exhibit it. The `skinnb=3.0` fix
   (v2.5.54) is therefore the correct *engine-level* remedy, not a sampling or
   barostat band-aid.
4. **Barostat divergence is deliberate**: FEP-SPell and FEP+ early stages use
   tighter coupling; this workflow moved to gentle Berendsen (`taup=5`) + larger
   skin to survive GPU NPT on a large complex box.

---

## 6. BAT.py (grounded from code)

Files read: `BAT/BAT.py`, `BAT/example-input-files/input-amber-short.in`,
`input-dd-rest-op.in`, `input-express-mbar.in`, `BAT/amber_files/*.in`,
`BAT/lib/{dd,dd-ti,sdr,sdr-ti,rest,equil}.py`.

BAT.py (Binding Affinity Tool) is an Amber/OpenMM ABFE/RBFE automation tool. It
is the most *route-flexible* of the four references: a single `fe_type` switch
selects the thermodynamic path.

### 6.1 Free-energy routes (`fe_type`, grounded)
- `express` — fast single-box double-decoupling (ranking; `input-amber-short.in`,
  `input-express-mbar.in`)
- `dd-rest` — double-decoupling with restraints (`input-dd-rest-op.in`)
- also `dd` (double decoupling), `sdr` (simultaneous decoupling-recoupling),
  `rest` (restraint-only), and relative (RBFE) variants in `lib/`
  (`dd.py`, `dd-ti.py`, `sdr.py`, `sdr-ti.py`, `rest.py`).

### 6.2 Alchemical components (grounded)
BAT.py decomposes ΔG into lettered **components**, each its own window set:
`m`/`n` (restraint attach in bound/free), `c` (conformational), `e`
(electrostatics/charge), `v` (van der Waals), `r` (release). Per-component step
counts (`input-amber-short.in`): equilibrium/production = m 25k/50k, c 25k/50k,
e 50k/50k, v 100k/100k steps.

### 6.3 λ / weight schedules (grounded)
- **Restraint attach/release** `attach_rest` = **10 windows**:
  `[0.00 0.10 0.24 0.56 1.33 3.16 7.50 17.78 42.17 100.0]` (force-constant
  weights, %).
- **Equilibrium gradual release** `release_eq` = **4**: `[10.0 2.5 0.5 0.0]`.
- **Decoupling λ** (express-MBAR) `lambdas` = **23 windows**:
  `[0.0001 0.02 0.04 0.06 0.08 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70
  0.75 0.80 0.85 0.90 0.92 0.94 0.96 0.98 0.9999]`.
- Integration: `dec_int = ti` (Gaussian-quadrature TI) **or** `mbar` —
  user-selectable, with `blocks = 5` for block-average error bars.

### 6.4 MD control (grounded)
`Temperature = 298.15`, `cut = 9.0`, `gamma_ln = 1.0`, **`barostat = 2`
(Monte-Carlo)**, **`dt = 0.004`** with **HMR `hmr = yes`**, eq steps
100k (gradual release) + 1,000k (post-release). Restraint force constants:
`rec_dihcf=50`, `rec_discf=5`, `lig_distance=5`, `lig_angle=250`,
`lig_dihcf=30`, `rec_com=lig_com=10` (kcal/mol per Å² or rad² as noted).

### 6.5 Relationship to this workflow
- **Restraint scheme differs**: BAT.py uses a **multi-component, multi-DOF**
  restraint set (distance + angle + dihedral + COM, attached over 10 weight
  windows) rather than this workflow's single Boresch `DISANG` introduced as one
  TI leg. BAT.py's is closer to APR-style staged restraint handling.
- **Barostat**: BAT.py defaults to **MC (`barostat=2`)**, like the reference platform/FEP+ — i.e.
  BAT.py would also *not* hit the pmemd.cuda fixed-grid "box changed too much"
  halt in the same way, because MC barostat volume moves + Amber's pairlist are
  handled differently than the Berendsen-NPT path this workflow used. (BAT.py
  also runs HMR `dt=0.004`, a different integration regime.)
- **HMR by default**: BAT.py ships `hmr=yes, dt=0.004`; this workflow uses
  `dt=0.002` without HMR. Adopting HMR is a possible future throughput win.
- **Route flexibility**: BAT.py's `express`/`dd`/`sdr`/`dd-rest` switch is more
  general than this workflow's fixed decharge→vdw→restraint leg structure;
  `sdr` (simultaneous decoupling-recoupling) is notable for avoiding a separate
  solvent leg.

### 6.6 Caveat
BAT.py bundles its **own vendored `pymbar`** (`lib/pymbar/`) — analysis numbers
from BAT.py come from that copy, which may differ in version from the `pymbar`
used elsewhere. Not a concern for protocol comparison, but relevant if
cross-checking ΔG estimator output.