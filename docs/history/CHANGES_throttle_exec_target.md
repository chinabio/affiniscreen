# Execution target + optional submission throttle

## Execution target (all four launch paths)
- Setup & Launch now has an **Execution target** selector: **GPU queue (LSF)**
  (default) or **Local host**. Applies to MM-GBSA (Amber + OpenMM), RBFE, ABFE.
- OpenMM MM-GBSA gained CLI flags: `--submit {local,gpu} --queue --walltime
  --project --conda-env --n-gpu`. `--submit gpu` writes a #BSUB script and
  bsubs it; `--submit local` (default) is the original behaviour.
- When GPU queue is selected, the queue/walltime are routed into the existing
  queue-aware paths (run_amber, fep_driver, OpenFE) without changing Local-host
  behaviour.

## Optional submission throttle (default OFF)
- New stdlib-only `amber_md/throttle_submit.py`: a detached login-node helper
  launched under nohup. Keeps no more than N of THIS batch's jobs in the queue
  (PEND+RUN), submitting the next as slots free. **Survives the page closing**
  (the previous GUI-side throttle did not).
- Exposed on Setup & Launch (MM-GBSA fan-out) and the FEP Campaign run tab
  (RBFE/ABFE) as a checkbox (default OFF) + "Max jobs in queue" number
  (default 8).
- DEFAULT BEHAVIOUR UNCHANGED: with throttle OFF, all jobs submit at once
  exactly as before.

## Verified
- submit_campaign default path: 12/12 submitted directly, no behavioural change.
- submit_campaign throttled path: 0 direct bsub, 12-job list handed to a
  DETACHED (start_new_session) throttler with --max-inflight honoured.
- throttle_submit stress test: peak concurrent in-flight == cap, never exceeded.
