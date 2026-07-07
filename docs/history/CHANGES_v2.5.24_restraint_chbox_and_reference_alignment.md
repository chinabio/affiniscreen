# v2.5.24 — ABFE restraint leg: align dual-copy build with validated FEP-SPell-ABFE

Follow-up to v2.5.23 (Option B dual-copy). The dual-copy *idea* was correct, but a
line-by-line comparison against the published, validated FEP-SPell-ABFE protocol
(`abfe/abfe_alchemy_morph.py`, `abfe/abfe_alchemy_md.py`,
`abfe/utils/common_tools.py::parse_leapin`, `abfe/md/md_tools.py::gen_masks`)
found 5 divergences. Two are likely fatal and explain the repeated
`complex_restraint` failures.

## Issues found (priority order)

### 1. FATAL (octahedral box): missing `ChBox` after the dual-copy tleap build
`SystemConfig.box_shape` defaults to **`octahedral`**, so `build_restraint_topology()`
calls `solvateOct`. The reference has an explicit, heavily-commented warning:

> When using solvateOct, `set Unit box {x,y,z}` only sets the box *dimensions*;
> the box alpha/beta/gamma angles get reset to 90 degrees, which is WRONG for a
> truncated octahedron and **causes failed simulations at the next step**.
> Fix: re-apply `ChBox` with the original `-al/-bt/-gm` angles (and X/Y/Z) on the rst7.

v2.5.23 never re-applied `ChBox`, so every octahedral restraint-leg run produced a
geometrically-wrong box → instability/crash on the very first MD step. **This is the
most likely cause of the four consecutive complex_restraint failures.**

**Fix:** after `saveAmberParm`, run `ChBox` on `complex_restraint.rst7`, restoring
`al/bt/gm` (=109.471° for a regular truncated octahedron) and the X/Y/Z box lengths
taken from the source equilibrated box.

### 2. Build from EQUILIBRATED coordinates, not the pristine build
The reference strips the protein from the **equilibrated** complex
(`equil_last_md_info['last_rst7']`) and re-uses the *equilibrated* box, only
rebuilding the box records. v2.5.23 strips from the pre-equilibration
`complex.prmtop`/`inpcrd` and runs a fresh `solvateOct` + fresh random ions, so the
restraint leg starts from an un-equilibrated box with new ion placement and must
re-equilibrate while the fixed Boresch restraint is already on → early blow-ups.

**Fix:** accept an optional `equil_inpcrd` (the complex `eq.rst`/last restart) and
strip the protein pose from THAT; carry the original box X/Y/Z into `ChBox`.

### 3. `crgmask=':2'` is an addition beyond the reference
The reference restraint leg is pure `:1`/`:2` TI with `ifsc=0` and **no crgmask** —
the two ligand copies are parameter-identical and only the Boresch potential is
perturbed. `crgmask=':2'` additionally perturbs electrostatics, deviating from the
validated cycle and risking double-counting vs. the decharge leg.

**Fix:** default `FEPConfig.restraint_crgmask = ""` (omit the crgmask line). Keep the
field so it can be re-enabled deliberately.

### 4. Do NOT re-run addIons on the dual-copy combine
With issue 2 fixed (inherit the equilibrated, already-neutral box) there is no need
to re-add ions, and re-adding them interacts badly with crgmask (over-neutralizing a
discharged `:2`). The reference does not re-run addIons in the restraint morph step.

**Fix:** default to NOT adding ions (`reion=False`). For a charged ligand built from
scratch, the guarded `ligand_charge` path remains available.

### 5. Boresch remap shift convention (verify, do not assume)
The reference shifts RECEPTOR masks by +1 *residue* (`gen_masks(shift=1)`) and leaves
ligand masks unchanged, because in THEIR complex the ligand is residue `:1`. In OUR
complex the ligand is the LAST residue, so `remap_boresch_for_dualcopy()` correctly
does the inverse atom-index mapping. It is internally consistent but cannot be
cross-checked against the reference numbers — validate against the actual built
`complex_restraint.parm7` with `tools/restraint_dualcopy_smoketest.py`.

## Files changed
| File | Change |
|---|---|
| `amber_md/abfe_restraint_topology.py` | Add `_chbox_octahedron()` + call it after tleap; accept `equil_inpcrd`; default `reion=False`; strip from equilibrated coords when given. |
| `amber_md/config.py` | `restraint_crgmask: str = ""` (was `":2"`). New `restraint_reion: bool = False`. |
| `amber_md/fep.py` | `_crgmask_block`: emit nothing for the restraint stage when `restraint_crgmask` is empty. |
| `amber_md/fep_driver.py` | Pass the equilibrated complex restart into `build_restraint_topology(..., equil_inpcrd=...)`. |

## Validation done (no Amber required)
- `box_shape` confirmed `octahedral` by default → ChBox gap was active on every run.
- Reference `parse_leapin(duplicate=True)` → `combine { Comp0 Comp0 protein }`,
  matching our `combine { LIG_R LIG_D protein }` order.
- `crgmask` / soft-core absent from all reference `COMPLEX_RESTRAINT_*` templates.

## Still requires the cluster
- tleap + ChBox build of `complex_restraint.parm7/.rst7`.
- pmemd.cuda smoke (λ=0.5): finite `VDWAALS`, finite `EPtot`, sane `DV/DL`, and a
  box that survives the first NPT step (the ChBox fix).
