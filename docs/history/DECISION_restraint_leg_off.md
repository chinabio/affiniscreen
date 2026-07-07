# Decision Record — Restraint MD leg OFF by default (v2.5.68 → v2.5.69)

**Status:** ADOPTED · **Date:** 2026-06-24 · **Applies to:** ABFE (Option-A) complex branch

---

## 1. TL;DR

The explicit ~2 ns Boresch **restraint MD leg** (`complex_restraint`) is **OFF by
default** starting v2.5.68. Its sampled free energy is **~0 by construction**, so
we skip the MD and add the **analytic** Boresch standard-state correction
(e.g. `-11.447503` kcal/mol) directly onto `complex_vdw`.

- `dG_bind` is **mathematically identical** with or without the leg (verified to
  machine precision: 5.447503 vs 5.447503 on the CI fixture).
- The Boresch **potential** is still **held ON** during `complex_decharge` and
  `complex_vdw` — the decoupled ligand still cannot drift. (CI-guarded invariant.)
- `--restraint-leg` restores the explicit MD leg for anyone who wants the
  empirical "ΔG_apply ≈ 0" validation (BAT.py / attach-pull-release style).

---

## 2. What the −11.447503 actually is (and where it comes from)

It is **not** measured from MD. It is a closed-form value written to
`boresch_correction.txt` by `amber_md/boresch.py::boresch_correction()`
(Boresch et al., *J. Phys. Chem. B* 2003, Eq. 32):

```
dG_release = kT * ln[ 8 π² V0 · sqrt(Kr·KθA·KθB·KφA·KφB·KφC)
                      ----------------------------------------- ]
                      ( r0² · sin θA · sin θB · (2π kT)³ )

dG_correction = − dG_release      (ADD direction; what we add to the cycle)
```

with `kT = 0.59219` kcal/mol @ 298 K, `V0 = 1660.54 Å³` (1 M standard state),
and Amber `rk → K = 2·rk`. We reproduced `-11.447503` to six decimals from the
package's own function using this ligand's measured geometry (`r0 ≈ 5.17 Å,
θ ≈ 100°`) with the standard force constants (`kr=10, kθ=kφ=100`). It depends
ONLY on the restraint definition + temperature, which is why it is analytic and
has **zero statistical error**.

---

## 3. What the restraint leg's MD does — and does NOT do

| Question | Answer (verified in source) |
|---|---|
| Does its MD contribute free energy? | **~0.** `dG_uncorrected_kcal_mol = 0.0`; total leg dG = the analytic term only. In Option-A the restraint is identical in every λ window (V0==V1), so the alchemical work integrates to ~0 by construction. |
| Does it provide starting coordinates (`prod.rst`) to decharge/vdw? | **No.** Confirmed: there is NO coordinate handoff anywhere in the code, and no dead/commented code that intended one. `complex_decharge` and `complex_vdw` each start from `build/complex.inpcrd` and run their **own** `min→heat→dens→eq→prod` ladder. |
| Does it keep the ligand in the pocket during decoupling? | **No — the restraint POTENTIAL does that, not the leg.** `fep.py`: "decharge / vdw stages: FIXED full k (restraint HELD ON the whole time so the ligand stays in the pocket while it is decoupled)." This is independent of whether the MD leg exists. |
| Why was the leg added (v2.5.16)? | For **formal cycle closure** — to give the analytic standard-state term an explicit home. The prior approach wrote the correction on `complex_vdw` with no MD leg and was called "formally incomplete." v2.5.68 returns to writing it on `complex_vdw`, but now with full analyzer support and an equivalence proof. |

**Myth checked and rejected:** "the restraint leg was secretly meant to generate
a better `prod.rst` to stop `complex_vdw` from failing." There is no such wiring,
no dead code, and the actual vdw-stability mechanism (restraint potential ON
during vdw + 1 fs dt + box-drift recovery) is unrelated to the leg.

---

## 4. Literature / other-workflow basis

| Method / code | Restraint contribution in the bound state | Separate MD leg? |
|---|---|---|
| **Boresch et al. 2003** (original) | analytic from force constants (the formula we use) | **No** |
| **Gilson / BAT.py** (attach-pull-release) | short "attach" ramp (the *a* component), deliberately small (~0) | Yes (short, validation) |
| **BFEE2/3, Gumbart–Roux–Chipot** | per-DOF PMF; small fully-interacting terms often (semi-)analytic | Mixed |
| **FEP+ ABFE / OpenFE-style automated ABFE** | analytic standard-state correction in the complex | **No** |

**Consensus:** the analytic standard-state term is the *mandatory* part; the
sampled "attach restraint to the bound, interacting complex" MD is the *optional*
~0 part — kept by some codes as a short check, skipped by others (Boresch
original, FEP+/OpenFE). Our default (`--no-restraint-leg`) matches the latter; our
opt-in (`--restraint-leg`) matches the former.

---

## 5. Why this is safe (and the guardrail)

Removing the leg only removes a ~0 sampled term and relocates the analytic
correction. The thing that protects `complex_vdw` (the restraint **potential**
during decoupling) is untouched. To make sure that can never be accidentally
broken, a CI invariant in `tools/check_dt_regression.py` asserts:

> "Boresch potential held ON during decharge/vdw" must remain present in
> `amber_md/fep.py`.

If a future edit removes the held-on restraint during vdw, **CI fails loudly.**

---

## 6. How to choose the mode

```
# DEFAULT (recommended): no restraint MD leg; analytic term folded onto complex_vdw
python -m amber_md.fep_driver ... --mode abfe          # --no-restraint-leg is implicit

# OPT-IN validation: run the explicit ~2 ns restraint MD leg as well
python -m amber_md.fep_driver ... --mode abfe --restraint-leg
```

Analysis is automatic in both layouts:

```
python tools/analyze_campaign.py <campaign_root> --recurse
```

It reads the analytic term from `complex_vdw/boresch_correction.txt` when there
is no restraint leg, and from `complex_restraint/` when there is — never both,
no double-counting.

---

## 7. Verification (CI fixture, v2.5.68)

```
dG_bind  with leg = 5.447503
dG_bind  no  leg  = 5.447503     identical = True ; both TRUSTED
CI: 24/24 PASS  (incl. 4 new guards: flags present, OFF-by-default,
                 vdw carries correction, restraint-potential-ON invariant)
```
