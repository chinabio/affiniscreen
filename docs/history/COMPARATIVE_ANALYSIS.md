# ABFE Workflow — Comparative Analysis & Design Rationale

**Document version:** 1.0 (paired with amber_md_workflow v2.5.55)
**Date:** 2026-06-21
**Scope:** How this Amber/pmemd ABFE workflow relates to the three external
references in the project — Schrödinger **FEP+** (the `.msj` protocol actually
run on this system), **the reference platform** (OpenMM ABFE), and **FEP-SPell-ABFE / BAT.py**
(the open-source Amber ABFE workflows this code is modeled on). All numbers
below are extracted from the actual artifacts in this project, not from
literature defaults. Where a value is *not* grounded in an uploaded file it is
explicitly marked as such.

---

## 1. Executive summary

| Workflow | Engine | Barostat | λ windows (complex) | Prod / window | Box-rebuild ("changed too much") failure mode? |
|---|---|---|---|---|---|
| **This workflow (Amber)** | pmemd.cuda | Berendsen (`barostat=1`, taup=5) | 21 (restraint leg) | 5 ns (2.5M × 2 fs) | **Yes** — fixed nonbond grid; cured in v2.5.54 via `skinnb=3.0` |
| **FEP+** (Schrödinger) | Desmond `desmond:fep` | MTK / NPT→μVT | **108** | **10 ns** (10000 ps) | No — Desmond re-tiles each step |
| **the reference platform** | OpenMM | Monte-Carlo barostat | (see the reference platform paper) | (see the reference platform paper) | No — OpenMM auto-regrids |
| **FEP-SPell / BAT.py** | pmemd / sander | MC / Berendsen | schedule per paper | per paper | Same engine family as this workflow |

The central engineering finding of this project (v2.5.52 → v2.5.54): the
"Periodic box dimensions have changed too much" halts that plagued the Amber
restraint leg are a **pmemd.cuda nonbond cell-list/skin artifact, not physical
box drift**, and not a failure mode the Desmond/OpenMM-based tools can even
exhibit. This is why FEP+ and the reference platform "just work" on the same chemistry while the
Amber port needed an explicit `skinnb` fix.

---

## 2. FEP+ protocol (extracted from the project `.msj` files)

Source files (all in this project):
`ATO_SEPTTR_absolute_binding_1.msj` (top-level workflow),
`..._complex.msj`, `..._solvent.msj`, `..._md.msj`,
`..._chg_complex.msj`, `..._chg_solvent.msj`, `..._chg_md.msj`,
and `ATO_SEPTTR_absolute_binding_1_multisim.log`.

### 2.1 Stage graph (from the top-level `.msj`)
1. `forcefield_builder_launcher` — OPLS parameterization
2. `fep_absolute_binding_md_launcher` — MD relaxation (`_md.msj`)
3. `fep_absolute_binding_fep_launcher` — the two FEP legs:
   - **complex** leg → `_complex.msj`
   - **solvent** leg → `_solvent.msj`
   - a separate **charge** pathway (`_chg_*.msj`) for the decharge step
4. `fep_absolute_binding_analysis` — cycle closure on the `.fmp` graph

### 2.2 λ schedule (grounded)
- **Complex leg: `lambda = "default:108"` → 108 windows** (`_complex.msj`, N = 108)
- **Solvent leg: `lambda = "default:60"` → 60 windows** (`_solvent.msj`, N = 60)
- Method: `task = "desmond:fep"`, `type = absolute_binding`

### 2.3 Production simulation (grounded, `_complex.msj` `lambda_hopping`)
- `time = 10000.0` ps → **10 ns production per window**
- `timestep = [0.004 0.004 0.008]` ps (RESPA; 4 fs bonded / 8 fs long-range)
- `ensemble = NVT` with **GCMC** water sampling
  (`gcmc.solvent_density = 0.03248`, `moves_per_cycle = 34000`,
  `mu_excess = -6.137`) → effectively μVT for buried-water exchange
- **REST / solute tempering** via `lambda_hopping` +
  `solute_tempering` (`exchange_probability = 0.3`, hot region by ASL),
  i.e. Hamiltonian replica exchange across the 108 windows
- `scale_solvent_vdw = 0.75` (the same 0.75 vdW scaling this Amber workflow
  uses — see `solvated_vdw` legs)

### 2.4 Pre-FEP relaxation (`_md.msj`, `task = "desmond:regular"`)
Brief staged relaxation ladder before FEP: `time = 100`, `12`, `20`, `20` ps
stages at `temperature = 10.0 → 298`, ending with a `time = 1000.0` ps (1 ns)
NPT equilibration. (Contrast with this workflow's min → heat → dens(0.5 ns) →
eq(1 ns) ladder per window.)

### 2.5 Observed runtime (grounded, `multisim.log`)
- Stage 3 (FEP) CPU time: **6h 04' 44"**; Stage 4 analysis: 1h 40' 51"
- **Total multisim duration: 12h 23' 56"** ("Multisim completed.")
- This is wall time on Schrödinger's GPU cluster for the full 108+60 window
  absolute-binding calculation of one ligand (lig 12944901 / ATO_SEPTTR).

> NOTE: the FEP+ "CPU time" lines are scheduler aggregate figures; the
> 12h 24m total is the figure to quote for end-to-end throughput.

---

## 3. This workflow (Amber) — what the project data shows

- Engine: **pmemd.cuda** (GPU), CPU `pmemd.MPI` fallback for density-settle.
- Restraint leg λ schedule (from the 21-window `prod.out` triage):
  `0.000, 0.004, 0.016, 0.036, 0.064, 0.100, 0.144, 0.196, 0.256, 0.324,
  0.400, 0.500, 0.600, 0.700, 0.800, 0.875, 0.925, 0.950, 0.975, 1.000`
  (21 windows) — i.e. **far coarser than FEP+'s 108**, by design (TI on a
  small, smooth restraint-introduction term needs fewer windows than a full
  REST-coupled decoupling).
- Production: `nstlim = 2,500,000` × `dt = 2 fs` = **5 ns/window**.
- Equilibration (grounded, `eq.out`/`dens.out`): eq 500k steps + dens 250k
  steps, completing identically for *every* window (ρ ≈ 1.011–1.012 after eq,
  ≈ 0.997 after dens) — including the 8 that later failed production.

---

## 4. The box-rebuild failure: root-cause analysis (the core of this project)

### 4.1 Symptom
8 of 21 restraint-leg windows (`0.064, 0.100, 0.144, 0.196, 0.256, 0.324,
0.400, 0.600`) died with the workflow logging "box drift" and `STAGE FAILED`
after a 10× regrid loop (including a CPU density-settle that itself reported
SUCCESS, then re-failed).

### 4.2 What the data actually showed (full 21-window triage)
- The failing `prod.out` files end on a **clean MD step** with **no error text
  in `prod.out` at all** — the "changed too much" message appears **only in
  `prod.console`** (pmemd stdout).
- At the moment of every halt the system is **physically healthy**: PRESS
  bounded (±90 bar), VOLUME varying **< 0.2 %**, density steady ≈ 1.016,
  **no NaN**, energies flat.
- Halts occur at wildly different step counts (0.1 %–22 % of nstlim) — there is
  **no consistent volume threshold**, which a genuine drift would show.
- `&cntrl` blocks are **byte-identical** between passing and failing windows;
  only the per-window Boresch `DISANG` differs.
- eq + dens end-states are **identical** for passing and failing windows →
  equilibration is not the cause.

### 4.3 Conclusion
This is the documented **pmemd.cuda nonbond cell-list / Verlet-skin error**, not
box drift. `skinnb` was unset (pmemd default **2.0 Å**) with `cut = 10.0 Å`; on
a large complex box the GPU grid rebuild intermittently fails on a *benign* NPT
fluctuation. The regrid recovery mdin **also lacked `skinnb`**, so every
"recovery" relaunched into the same lottery → the 10× loop and false failure.

```
benign NPT fluctuation
        │
        ▼
GPU nonbond cell-list rebuild lands wrong (skinnb too small)
        │
        ├─► prod.out: clean, truncated mid-run (no error inside)
        └─► prod.console: "Periodic box dimensions have changed too much"
                │
                ▼
        workflow greps console → logs "box drift" → regrid from eq.rst
                │  (regrid mdin ALSO lacks skinnb)
                └────────────► same failure → 10× loop → STAGE FAILED
```

### 4.4 Why FEP+ and the reference platform never see this
- **FEP+ / Desmond** re-tiles its spatial decomposition every step and does not
  use Amber's fixed pairlist-grid box check; the error class does not exist.
- **the reference platform / OpenMM** uses a Monte-Carlo barostat and auto-regrids the nonbonded
  grid on volume change; again, no fixed-grid halt.
- Only the **Amber/pmemd** family exposes this, which is exactly why the
  open-source Amber references (BAT.py, FEP-SPell) and this workflow must set an
  adequate `skinnb` — and why the earlier "make the barostat gentler" theory,
  while harmless, was not the real cure.

### 4.5 Fix shipped (v2.5.54)
`skinnb = 3.0` added to **prod (TI + restraint), eq, dens, AND the regrid
recovery mdin**. The v2.5.52 `taup = 5.0, barostat = 1` change is retained as a
complementary, benign stabilizer (it matches the eq/dens protocol).

---

## 5. Design comparison table

| Dimension | This workflow (Amber) | FEP+ (grounded `.msj`) | the reference platform (OpenMM) |
|---|---|---|---|
| Sampling enhancement | independent windows (TI/MBAR) | **REST + λ-hopping** (HREX), 108 win | per the reference platform paper |
| Water sampling | standard TIP3P NPT | **GCMC** buried-water exchange (μVT) | per the reference platform paper |
| Barostat | Berendsen (taup=5) | MTK NPT → NVT prod | Monte-Carlo |
| Prod length | 5 ns/window | **10 ns/window** | per the reference platform paper |
| Restraint | Boresch (DISANG), TI-introduced | absolute_binding restraints, fc=50 | Boresch-type |
| vdW solvent scaling | 0.75 | **0.75** (`scale_solvent_vdw`) | per the reference platform paper |
| End-to-end wall time | (per the cluster run) | **~12h 24m** (1 ligand) | per the reference platform paper |
| Fixed-grid box halt | **possible → fixed (skinnb)** | n/a | n/a |

---

## 6. Caveats / honesty notes

1. **the reference platform numbers are intentionally left as "per paper"** here. The the reference platform code
   (`reference-platform.zip`) and paper (`2603.22274v1.pdf`) are in the project, but I
   have not re-extracted a verified λ-count / per-window time table for the reference platform in
   this document; do not quote specific the reference platform windows from this file.
2. The FEP+ figures **are** grounded in the uploaded `.msj`/`multisim.log` and
   can be cited directly.
3. The 12h 24m FEP+ total is for **one ligand** on Schrödinger hardware; it is
   not directly comparable to a per-window Amber cost without normalizing by
   window count and GPU type.
4. The `skinnb = 3.0` value is the conservative, widely-used remedy; if a window
   still trips, escalate to 3.5–4.0 (documented in the v2.5.54 changelog).
