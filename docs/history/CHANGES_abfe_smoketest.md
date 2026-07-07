# v2.5.0 build final32: fast single-ligand ABFE smoke-test launcher

## New: run_abfe_smoketest.sh
ONE small ABFE job sized to finish fast, to validate the whole pipeline on
the login node/the cluster before spending campaign hours. RUN ONLY AFTER preflight GO.

  cd ~
  bash $WF/run_abfe_smoketest.sh protein.pdb ligand.sdf [RESNAME]
  bjobs
  cat <WORKDIR>/fep/ABFE_RESULT.txt

Design (matches fep_driver v2.4.21 real CLI):
  * --mode abfe with --protein-pdb/--ligand-file -> auto topology build.
  * --auto-boresch -> Boresch restraints from the built complex PDB.
  * Driver runs on the LOGIN node and --submit bsubs each lambda window +
    per-leg analyzers + a cycle-closer that writes fep/ABFE_RESULT. The
    launcher is NOT a gpu -Is job (the gpu queue refuses interactive jobs).
  * cd $HOME first -> avoids the broken-CWD 'module' error.
  * --project your-project explicit (also the default now).

FAST sizing (too short for a real dG; goal = pipeline closes):
  --nstlim-eq 10000 (20 ps), --nstlim-prod 50000 (100 ps),
  5-point lambda schedules per stage, --walltime 2:00/window. Expect noise.

Resume after preemption:
  python -m amber_md.fep_driver --resume --work-dir <WORKDIR> --submit \
         --project your-project
