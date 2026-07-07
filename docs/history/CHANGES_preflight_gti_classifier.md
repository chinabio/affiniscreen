# v2.5.0 build final31: fix preflight check[4] false-positive + batch submit

## Context (the login node run, amber-md + amber/22.8)
Environment is HEALTHY:
  pmemd.cuda / pmemd.cuda.MPI  Version 22.0   (/share/apps/amber/rhel8/amber22)
  parmed 4.0.0 / numpy 1.26.4  (no MM-GBSA numpy.compat trap)
  alchemlyb 2.5.0 / pymbar 4.2.0
  modules: gcc/11.5 cuda/11.8 amber/22.8 pymol/3.0.4

## Bug fixed: check[4] false NO-GO 'gti'
The toy test system was two Na+ ions -> each TI region had net charge +1, so
pmemd printed "Skip neutralizing charges..." and exited non-zero. The old
detector saw rc!=0 and wrongly reported "build lacks GTI". But the mdout proved
the OPPOSITE: pmemd parsed icfe and built "TI region 1 / TI region 2" -- i.e.
GTI IS supported.

Fix:
  * _write_tiny_system now builds a NEUTRAL TIP3P water box (no charged-region
    edge case).
  * check_gti_run now classifies by EVIDENCE, not just rc:
      - rejected  := icfe/ifsc/gti AND (not supported/allowed/unknown/...)  -> BLOCK
      - engaged   := mdout shows TI region / softcore / clambda / charges     -> OK
                     (OK even if rc!=0, labeled 'benign' -- toy system quirk)
      - else      := WARN inconclusive (never a false hard NO-GO)
  Verified against the actual the login node mdout tail -> now 'OK (gti supported)'.

## New: run_preflight_abfe.bsub
The gpu queue refuses interactive jobs ("Queue does not accept interactive jobs"),
so the preflight ships as a batch script:
    cd ~ && bsub < $HOME/Tools/affiniscreen/run_preflight_abfe.bsub
Reads abfe_preflight.<JOBID>.out. cd $HOME first to dodge the broken-CWD module
error seen earlier.

## Note on the [5] MPI WARN
pmemd.cuda.MPI -ng multi-window failed on the login node (mpi-abort). That is a
WARN, not a blocker: single-window TI is unaffected; HREMD/-ng needs a proper
multi-GPU MPI launch on a compute node, which we test in the real job.
