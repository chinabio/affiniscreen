# amber_md v2.5.69 — clean-run release (restraint leg OFF by default)

This is a **documentation + clean-run-readiness** release on top of v2.5.68.
No functional change to the free-energy math; ships the decision record, an
updated clean-run guide, and version-file sync so a fresh ABFE campaign can be
launched and analyzed end-to-end.

## Carried from v2.5.68 (functional)
- **Restraint MD leg OFF by default.** `complex_restraint` is no longer generated.
  The analytic Boresch standard-state correction is folded onto `complex_vdw`.
- **`--restraint-leg`** restores the explicit ~2 ns MD leg (BAT.py-style check).
- **INVARIANT (CI-guarded):** Boresch potential held ON during decharge+vdw, so
  the decoupled ligand cannot drift in either mode.
- `analyze_campaign.py` supports both layouts; `dG_bind` identical to machine
  precision (verified 5.447503 == 5.447503; 24/24 CI PASS).

## New in v2.5.69 (docs / housekeeping)
- **`DECISION_restraint_leg_off.md`** — full rationale: what −11.4475 is, what the
  leg does/does not do (with the "no coordinate handoff" and "potential held on"
  findings), literature basis (Boresch 2003, BAT.py, BFEE, FEP+/OpenFE), safety
  argument, and the CI invariant.
- **`README_CLEAN_RUN.md`** updated to v2.5.69: documents the default no-leg
  layout, the `--restraint-leg` opt-in, and the pre-submission checklist.
- **`VERSION`** file synced to match `amber_md/__init__.py` (was stale at 2.5.61).
- CI banner updated.

## Clean-run quick start
```
python -m amber_md.fep_driver ... --mode abfe          # no restraint leg (default)
python tools/analyze_campaign.py <campaign_root> --recurse
```

## Next (after your clean run)
Collect `complex_vdw/lambda_*/prod.out`, `prod.in`, and any `prod.console.*` /
LSF `.err` from a failing window so the `complex_vdw` instability (soft-core
endpoint overflow / dt / box-drift) can be diagnosed on real data. This is
orthogonal to the restraint-leg change.
