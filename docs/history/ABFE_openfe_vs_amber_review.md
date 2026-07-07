# Improving the Amber ABFE path using the (working) OpenFE ABFE as reference

Scope: compared `amber_md/abfe_openfe_plan.py` (OpenFE `AbsoluteBindingProtocol`,
known-working) against the Amber/pmemd ABFE path
(`amber_md/fep.py` + `amber_md/fep_driver.py` + `amber_md/boresch.py` +
`amber_md/config.py`). Goal: find where the Amber path can be brought in line
with the validated OpenFE protocol.

## Side-by-side

| ABFE element | OpenFE (working) | Amber path | Gap |
|---|---|---|---|
| Decoupling stages | charge annihilation + LJ decoupling | `decharge` + `vdw` two-stage | aligned |
| Soft-core | internal, validated | GTI `gti_add_sc=1`, `scalpha=0.5`, `scbeta=12` | aligned |
| Boresch restraint | **ramped via a dedicated `lambda_restraints` leg**; MBAR measures the attach work | **full-strength restraint on EVERY complex window**; single analytical RRHO term | **BIAS — #1** |
| Sampler | **HREMD (replica exchange) by default** | independent windows by default; HREMD optional | **#2** |
| Restraint atom pick | receptor anchor gated to `host_min/max` (0.5–1.5 nm) to avoid PBC | nearest-Cα by **index arithmetic** (`a+3`, `a+6`) | **#3** |
| Box stability (eq) | MC barostat, validated | MC barostat caused λ=0.50 blow-up — fixed in final60 (gentle eq + temp gate) | addressed |

The two-stage decharge/vdw decomposition and the soft-core settings are sound.
The substantive gaps are **restraint handling** and **sampling** — exactly the
parts OpenFE's docstring says it "does internally and is heavily validated."

---

## #1 (HIGH, correctness) — restraint attach work is unaccounted → systematic dG_bind bias

**Evidence in code.**
- `fep.py::setup_leg` loops over *all* `_active_lambdas` and calls
  `_write_boresch_RST(wd)` whenever `self.boresch` is set (lines 422–429).
- `fep.py::_write_boresch_RST` emits the same `boresch.RST` with **constant**
  `rk2/rk3` (= `b["kr"], b["kth"], b["kph"]`) — no λ scaling (lines 390–408).
- `boresch.py::boresch_correction` is the standard RRHO expression for
  **removing a fully-on restraint from the fully-decoupled end state**.

**Why it biases the result.** The full restraint is present while the ligand is
still (partly) coupled — throughout `complex_decharge` and the coupled end of
`complex_vdw`. The free-energy cost of **attaching** the restraint to the
*interacting* ligand is neither sampled nor included in the analytical term
(which is exact only for the non-interacting, fully-restrained end state).
OpenFE measures this with its `lambda_restraints` leg, so its dG_bind is
complete. The Amber number will differ from the validated OpenFE value by the
unaccounted attach work — even when every leg "completes" cleanly. This is the
most likely cause of an Amber/OpenFE ABFE discrepancy.

**Proposed fix (mirror OpenFE).** Add a restraint-attach sub-leg to the complex
ABFE schedule:
- New leg `complex_restraint` with its own λ grid (`restraint_lambdas`,
  e.g. `(0.0,0.1,0.2,0.35,0.5,0.65,0.8,0.9,1.0)`).
- In that leg the ligand is **fully coupled** (no decharge/vdw decoupling) and
  the Boresch force constants are scaled by λ: `rk2=rk3=λ·k_full`. Implement by
  having `_write_boresch_RST` accept a `scale` argument driven by `clambda`.
- MBAR/TI over that leg yields `dG_attach`. The total becomes
  `dG_bind = -(dG_attach + dG_decharge^cplx + dG_vdw^cplx)
             + (dG_decharge^solv + dG_vdw^solv) + dG_restraint_RRHO`,
  where the analytical RRHO term remains for the **fully-decoupled** removal
  (its valid regime). `absolute_binding_dG_from_legs` must be extended to take
  the attach leg.
- This is a **protocol change** (extra leg + accounting), so it is gated behind a
  config flag `abfe_ramp_restraint: bool = True` and validated on the toluene/
  T4-lysozyme system OpenFE ships, target agreement <~0.5 kcal/mol.

> NOTE: this is a design change, not a one-line patch. Recommend implementing on
> a branch and validating numerically before making it the default.

## #2 (MEDIUM, convergence) — default the ABFE complex legs to HREMD

OpenFE's default sampler is replica exchange. The Amber default is independent
windows (`hremd=False`); HREMD is only opt-in. The soft-core overlap region —
where independent sampling under-converges and where λ=0.50 detonated — is
exactly where REX helps. Now that **both** eq paths (LSF-array and HREMD) carry
the final60 stability gate, enabling HREMD for ABFE complex legs is low-risk.

**Proposed change (low-risk):** in `fep_driver.py`, when `mode == "abfe"`,
default `a.hremd = True` for the two `complex_*` legs unless the user passed
`--no-hremd`. Solvent legs can stay independent (cheap, well-behaved). Add
`abfe_hremd: bool = True` to `FEPConfig`.

## #3 (LOW, robustness) — PBC-safe Boresch receptor anchor selection

`boresch.py::select_boresch_atoms` picks the receptor reference atoms by Cα
**index offset** (`a_idx+3`, `a_idx+6`), which can land on spatially distant or
near-box-edge atoms → unstable restraint geometry. OpenFE constrains the anchor
distance to `[host_min_distance, host_max_distance]` (0.5–1.5 nm).

**Proposed change:** choose `b` and `c` by *distance* from `a` within a
`[host_min_A, host_max_A]` shell (defaults 5–12 Å) and require a minimum
pairwise separation so the three define a stable frame, instead of fixed index
offsets. Add `host_min_A`, `host_max_A` to config. Pure setup-time change; no
effect on the MD protocol.

---

## Recommended order
1. **#3** first — small, self-contained, improves restraint stability immediately
   and reduces the chance of future λ blow-ups at the restraint anchor.
2. **#2** next — one-flag default change, leans on the final60 eq gate already in
   place; improves convergence to match OpenFE's sampler.
3. **#1** last and on a branch — the real correctness fix, but it is a protocol
   change that must be numerically validated against the OpenFE reference number
   before becoming default.

## What NOT to change
- Two-stage decharge/vdw decomposition and soft-core (`scalpha/scbeta`,
  `gti_add_sc`) — already aligned with OpenFE; leave as-is.
- The analytical RRHO term — keep it for the fully-decoupled removal; it is
  correct there. #1 *adds* the attach leg, it does not replace the RRHO term.
