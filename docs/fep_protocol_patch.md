# amber_md v2.5.16 — ABFE protocol correctness patch

Derived from a line-by-line comparison of `amber_md/fep.py` (v2.5.15) against
the **published, validated** FEP-SPell-ABFE protocol
(`abfe/md/amber_alchemy_mdin.py`, J. Chem. Inf. Model. 2025, 65, 2711).

This patch addresses why ABFE windows crash (box drift / eq blow-up / T≈15000 K
in the λ=0.6–0.9 soft-core "danger zone") and two systematic free-energy biases.

---

## 1. GTI soft-core control block (THE crash fix)

### Current `amber_md/config.py` (FEPConfig) — replace these defaults
```python
    scalpha: float = 0.5          # -> 0.2
    scbeta: float = 12.0          # -> 50.0
    gti_add_sc: int = 1           # keep
    gti_chg_keep: int = 0         # -> 1
```

### Add these new FEPConfig fields (full GTI soft-core set, BAT.py/FEP-SPell-proven)
```python
    # v2.5.16: complete GTI soft-core controls. Without the ele/vdw/lam_sch
    # trio the TI soft-core path is half-configured -> stiff dV/dl -> the
    # mid-lambda blow-ups. These match the published FEP-SPell protocol.
    scalpha: float = 0.2          # was 0.5  (AMBER-recommended soft pair)
    scbeta: float = 50.0          # was 12.0
    gti_lam_sch: int = 1
    gti_ele_sc: int = 1
    gti_vdw_sc: int = 1
    gti_scale_beta: int = 0
    gti_cut_sc: int = 0
    gti_cut: int = 1
    tishake: int = 1              # SHAKE on the TI region (stability)
    gti_chg_keep: int = 1         # was 0
    logdvdl: int = 0
```

### Replace `FEP._softcore_block()` in `amber_md/fep.py`
```python
    def _softcore_block(self):
        """Complete GTI soft-core + smoothing block (v2.5.16).

        Matches the published FEP-SPell-ABFE protocol. The previous block
        emitted only scalpha/scbeta + gti_add_sc/gti_chg_keep, leaving
        gti_ele_sc / gti_vdw_sc / gti_lam_sch unset -> the half-configured
        TI soft-core that diverged at mid-lambda.
        """
        c = self.cfg
        g = lambda name, d: getattr(c, name, d)
        return (
            f"  scalpha={g('scalpha',0.2)}, scbeta={g('scbeta',50.0)},\n"
            f"  gti_lam_sch={g('gti_lam_sch',1)}, gti_ele_sc={g('gti_ele_sc',1)}, "
            f"gti_vdw_sc={g('gti_vdw_sc',1)},\n"
            f"  gti_scale_beta={g('gti_scale_beta',0)}, gti_cut_sc={g('gti_cut_sc',0)}, "
            f"gti_cut={g('gti_cut',1)},\n"
            f"  gti_add_sc={g('gti_add_sc',1)}, gti_chg_keep={g('gti_chg_keep',1)}, "
            f"tishake={g('tishake',1)}, logdvdl={g('logdvdl',0)},\n"
        )
```

---

## 2. vdW λ schedule — densify the REAL danger zone (0.575–0.80)

The published vdW/decharge-of-complex schedule is **44 windows**, with the
finest spacing exactly where amber's windows die. Replace `vdw_lambdas`:

```python
    # v2.5.16: published FEP-SPell vdW decoupling schedule (44 windows).
    # Ultra-dense across 0.575-0.80 where the ~60-85%-decoupled core lets
    # atoms nearly overlap and dV/dl turns stiff. This is the empirically
    # stable grid; the old 28-window grid was densest in the wrong region.
    vdw_lambdas: tuple = (
        0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.475,
        0.50, 0.525, 0.55, 0.575, 0.585, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64,
        0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75,
        0.76, 0.77, 0.78, 0.79, 0.80, 0.825, 0.85, 0.875, 0.90, 0.95, 1.0)
```

Decharge stays at the published 11-point grid (charges removed with vdW core
fully present, so it does not need the danger-zone density):
```python
    decharge_lambdas: tuple = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                               0.6, 0.7, 0.8, 0.9, 1.0)
```

> **Why this matters numerically:** with Δλ ≤ 0.01 across 0.575–0.80, each
> window's perturbation ΔU is small, so (a) the integrator stays stable and
> (b) nearest-neighbour MBAR overlap rises above the 0.03 floor the analyzer
> warns about — i.e. it *also* fixes the "MBAR ill-posed / u_nk rank" failures.

---

## 3. Add the dedicated restraint-removal leg (formal ABFE cycle)

The current code carries the Boresch potential on both complex sub-stages and
only ever subtracts the analytic correction. The published cycle instead runs a
**separate restraint leg** that turns the Boresch restraint ON (λ=0) → OFF
(λ=1) on the *interacting* complex, then applies the analytic standard-state
term at the decoupled end. 16-window schedule:

```python
    restraint_lambdas: tuple = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75,
                                0.8, 0.85, 0.9, 0.925, 0.95, 0.96,
                                0.97, 0.98, 0.99, 1.0)
```

Wire it into the two-stage leg list in `fep_driver.py` alongside
`complex_decharge` / `complex_vdw`, using the SAME mdin templates but with
`nmropt=1 / DISANG=boresch.RST` and `icfe=1` ramping the restraint via
`clambda` (the restraint leg has no soft-core: `ifsc=0`). The cycle-closer then
adds dG_restraint(MD) + dG_restraint(analytic, this code).

---

## 4. Equilibration: hold positional restraints through min+heat+press

Match the published staged equilibration so the soft-core ligand cannot drift
before production:
- minimization: `ntr=1, restraint_wt=5.0, restraintmask="<solute heavy atoms>"`
- heat (NVT) + press (NPT): keep `ntr=1` with the same mask, Boresch `nmropt=1`
- production: `ntr=0` (Boresch still on via `nmropt=1`)

amber's `_eq_in`/`_dens_in` already use the gentle `dt=0.001` + Berendsen path
(good); the missing piece is the **positional restraint mask** carried from
min through press.

---

## 5. Drop-in analysis fixes (already validated in this thread)

- **`boresch.py`** — corrected `rk`→`2·rk` factor. Self-test now returns
  **11.629 kcal/mol** vs. the published reference **11.62** (was 10.40, a
  +1.2 kcal/mol bias on every ABFE). Run `python -m amber_md.boresch --self-test`.
- **`charge_correction.py`** — new Rocklin PB finite-size correction
  (NET_USV + RIP + DSC) for net-charged ligands; call it on the final frame of
  `complex_decharge` and `solvent_decharge` and add the difference to dG_bind.

---

## Validation checklist after applying

1. `python -m amber_md.boresch --self-test`  → `PASS` (11.629).
2. Generate one `prod.in` and confirm it contains `gti_ele_sc=1, gti_vdw_sc=1,
   gti_lam_sch=1, scalpha=0.2, scbeta=50.0, tishake=1`.
3. Count windows: decharge=11, vdw=44, restraint=16.
4. Smoke-test one mid-zone window (λ≈0.65 vdW) — it should no longer hit the
   eq temperature gate (`eq_temp_max_K`).
5. For a charged ligand, confirm `charge_correction.csv` is written and the
   Total term enters `ABFE_RESULT.json`.
