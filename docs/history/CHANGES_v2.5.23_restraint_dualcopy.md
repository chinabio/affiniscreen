# v2.5.23 — ABFE restraint leg: dual-copy TI topology (Option B)

## Summary
Fixes the four consecutive failures of the formal-ABFE `complex_restraint` leg
(v2.5.20–v2.5.22) by giving the leg a **dedicated dual-copy topology** and
running it as **ordinary TI** instead of trying to ramp a Boresch potential on
the shared `complex.prmtop`.

## Root cause
All of v2.5.20–v2.5.22 reused `build/complex.prmtop`, where the ligand is the
**last** residue (`:872`, after 871 protein residues). With both TI end states
sharing the same fully-interacting ligand, every mask/ifsc combination failed:

| ver | restraint approach | failure |
|---|---|---|
| 2.5.20 | empty TI masks | `timask1/2 must match same #atoms` |
| 2.5.21 | `:LIG/:LIG`, ifsc=0 | hard-core `VDWAALS=****` overflow |
| 2.5.22 | `:LIG/:LIG` + soft-core (ifsc=1) | soft-core `SC_VDW=****` overflow |

## Fix (Option B, FEP-SPell-ABFE dual-copy)
Build `build_restraint/complex_restraint.parm7/.rst7` laid out as
`combine { LIG LIG protein }`:
- `:1` = real ligand (TI end state V0)
- `:2` = dummy ligand, charges zeroed via `crgmask=':2'` (TI end state V1)
- protein shifted by one ligand copy (+`n_lig` atoms)

The leg then runs **ordinary TI**: `timask1=':1'`, `timask2=':2'`, `ifsc=0`,
plus a **fixed-k** (non-λ-scaled) Boresch restraint. Matching 37-atom TI regions
→ no atom-count abort; no soft-core → no SC_VDW overflow; full ligand vdW never
goes through the TI hard-core path → no VDWAALS overflow.

## Files changed
| File | Change |
|---|---|
| `amber_md/abfe_restraint_topology.py` | **NEW.** `build_restraint_topology()` (strip protein from complex → tleap dual-copy → fresh solvate/ions) and `remap_boresch_for_dualcopy()` (shift Boresch atom indices into the ligand-first layout). |
| `amber_md/fep.py` | `_mask_block` restraint stage → `(':1', ':2', '', '')`; `_ifsc_value` restraint → `0`; `_crgmask_block` restraint → `':2'`; `_write_boresch_RST` force constants no longer λ-scaled. |
| `amber_md/fep_driver.py` | Builds the dual-copy topology once (reused on `--resume`), remaps the Boresch atoms, and points **only** the `complex_restraint` leg at the new prmtop. Decharge/vdw legs unchanged. Graceful fallback to `complex.prmtop` if mol2/frcmod unknown. |
| `amber_md/config.py` | New `FEPConfig` fields: `restraint_timask1=':1'`, `restraint_timask2=':2'`, `restraint_crgmask=':2'`, `build_restraint_topology=True`. |

## Charge neutrality
A **neutral** ligand (verified: net charge ≈ 0 for the test ligand) adds zero net
charge when copied, so the existing 15 Na⁺ / 15 Cl⁻ recipe is unchanged. For a
**charged** ligand the builder neutralizes `2 × ligand_charge` (both copies) via
`SystemConfig.ligand_charge`. Confirm the tleap "unperturbed charge ~0.000" line.

## Validation done (no Amber required)
- `complex.prmtop` parse: LIG = `:872`, atoms 13358–13394 (n_lig=37), n_prot=13357.
- `count_residue_atoms()` → (37, 872) ✓
- `remap_boresch_for_dualcopy()` → ligand atoms map into `:1` (1/18/37), protein
  atoms shift by +37 ✓
- All four edited/new modules pass `ast.parse` ✓

## Still requires the cluster
- tleap build of `complex_restraint.parm7/.rst7`
- pmemd.cuda smoke test (λ=0.5: finite `VDWAALS`, sane `DV/DL`) then full leg
- Final charge-neutrality confirmation after the dual-copy combine

## Backward compatibility
`build_restraint_topology=True` by default. Set it `False` to restore the old
(broken) behaviour. Solvent/decharge/vdw legs and the analytic Boresch
correction path are unchanged; the analyzer reads `complex_restraint/ti*.out`
as ordinary TI.
