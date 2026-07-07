# CLEAN-RUN GUIDE (amber_md v2.5.77)

## ABFE complex branch layout (v2.5.68+)

By DEFAULT the restraint MD leg is **OFF**. The complex branch legs are:

  - `complex_decharge`   (restraint POTENTIAL held ON)
  - `complex_vdw`        (restraint POTENTIAL held ON; carries the analytic
                          Boresch standard-state correction in
                          `complex_vdw/boresch_correction.txt`)

The Boresch standard-state term (e.g. `-11.447503` kcal/mol) is **analytic** —
it is NOT sampled. See `DECISION_restraint_leg_off.md` for the full rationale.

Opt in to the explicit ~2 ns restraint MD leg (BAT.py-style validation) with
`--restraint-leg`; that re-adds `complex_restraint/` and moves the analytic
term there instead. `dG_bind` is identical either way.

Pipeline per window (every leg):
  min -> heat (NVT ramp, restrained, GATED) -> dens (NPT, restrained)
      -> eq (NPT, restrained, GATED) -> prod (NPT; Boresch potential ON for
         decharge/vdw, free otherwise)

## Launch

```
# DEFAULT (recommended): no restraint MD leg
python -m amber_md.fep_driver ... --mode abfe

# OPT-IN validation: also run the explicit restraint MD leg
python -m amber_md.fep_driver ... --mode abfe --restraint-leg
```

## Before submitting
1. Generate legs (fep_driver).
2. `bash tools/verify_stage_order.sh <leg>` for `complex_vdw`,
   `complex_decharge`, `solvent_vdw`, `solvent_decharge` -> expect PASS.
3. First cluster check: `complex_vdw` lambda=0.5 ->
   min->heat->dens->eq->prod all exit 0, peak T < gates, finite
   EPtot/VDWAALS, no `SC_VDW=****`.
4. (Only if `--restraint-leg`) `python tools/check_restraint_leg.py
   <run_dir>/complex_restraint`  (Option-A reliability check).

## Analyze (handles BOTH layouts automatically)
```
python tools/analyze_campaign.py <campaign_root> --recurse
```
Reads the analytic Boresch term from `complex_vdw/boresch_correction.txt` when
there is no restraint leg, or from `complex_restraint/` when there is. Never
both — no double-counting.

## Clean-run-safe defaults (FEPConfig)
do_heat=True, heat_nstlim=100000, heat_dt_ps=0.001, heat_T_start=5.0,
heat_ramp_frac=0.8, heat_restraint_wt=5.0, heat_temp_max_K=1000.0, vlimit=20.0,
posres_mask_default='!:WAT,ions & !@H=', eq_dt_ps=0.001, eq_barostat=1,
restraint_crgmask='', restraint_reion=False. prod dt = 0.001 ps (1 fs);
prod = 10 ns. Env overrides: HEAT_TEMP_MAX_K, EQ_TEMP_MAX_K.

## Known follow-up
If `complex_vdw` windows fail, collect from a failing window:
`prod.out` (or `mdout`), `prod.in`, and `prod.console.*` / LSF `.err`.
The AMBER error line (e.g. `vlimit exceeded`, SHAKE failure, `* NaN *`,
soft-core `SC_VDW=****`) pinpoints the cause. This is independent of the
restraint-leg setting.

Version: amber_md v2.5.77. Single source of truth:
amber_md/__init__.py __version__; mirrored in VERSION and run_amber.py.
