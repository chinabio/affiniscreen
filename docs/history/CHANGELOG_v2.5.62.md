# amber_md v2.5.62

## Fixed -- production timestep regression (real cause of the 2 fs ABFE crashes)
v2.5.61 set config dt_ps/prod_dt_ps = 0.001, but 2 fs still reached generated
prod.in (dt=0.002, nstlim=5000000), blowing up ABFE complex_restraint mid-band
windows (lambda ~0.10-0.40; lig_12944901 lambda=0.100/0.144).

1. fep_driver.py  : --dt default 0.002 -> 0.001 (it overrode FEPConfig.dt_ps);
                    --nstlim-prod 2_500_000 -> 5_000_000.
2. config.py      : prod_nsteps 5M -> 10M ; nstlim_prod 1M -> 2M.
3. gui/amber_config.py : ns->steps /2 -> /1 (prod & equil).
4. md_inputs.py   : heat/equil dt 0.002 -> 0.001.
5. equilibration_fepspell.py / abfe_integration.py : timestep 0.002 -> 0.001.
6. gui/Home.py   : release notes updated; prose 'dt=0.002' literal removed.

## Added / Enforced
- amber_md/_dt_guard.py : assert_prod_dt_safe(cfg) -> UnsafeTimestepError if any
  effective production dt > 0.001 ps.
- WIRED into fep_driver._build_configs (chokepoint for run_fep + run_rbfe).
- tools/check_dt_regression.py : CI self-check (exit!=0 on regression).

## Verify
    python -c 'import amber_md; print(amber_md.__version__)'   # 2.5.62
    python -m amber_md._dt_guard                               # [dt-guard] OK
    python tools/check_dt_regression.py                        # RESULT: PASS
