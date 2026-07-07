# v2.5.3 build final54: parallel MM/GBSA scoring (1 GPU + N CPU cores)

## Request
Keep ONE GPU for pmemd.cuda MD, but grant N CPU cores alongside it so the
in-job MM/GBSA step runs MMPBSA.py.MPI in parallel (the scoring loop is
embarrassingly parallel over frames).

## What was already in place
submit.py._mmgbsa_tail_block already scales correctly at RUNTIME:
  NP=${LSB_DJOB_NUMPROC}; (cap via mmgbsa_n_cpu); if NP>1 -> mpirun -np NP MMPBSA.py.MPI
So scoring parallelism just needed LSF to actually GRANT cores.

## The hard part (the cluster `gpu` queue: slots == GPUs)
Plain `-n 8` reserves 8 GPUs (the final53 bug). To get cores WITHOUT extra
GPUs, request the GPU as a RESOURCE and let -n mean cores.

## Change
config.py HPCConfig:
  * gpu_request_mode: str = "rusage"   (default) | "slots" (legacy)
  * gpu_cpu_cores: int = 8             cores requested for the GPU job
submit.py._header() GPU branch, three paths:
  * mode="rusage" (DEFAULT):
        #BSUB -n {gpu_cpu_cores}
        #BSUB -R "rusage[mem=...,ngpus_physical={n_gpu}]"
        #BSUB -gpu "num={n_gpu}:j_exclusive=no"
    -> cores via -n, GPU via resource; cores DO NOT inflate GPU count.
  * mode="slots" (legacy, slots==cores clusters):
        #BSUB -n {n_gpu_slots} ; -gpu "num={n_gpu}"   (final53 behavior)
  * explicit n_slots arg (FEP windows): honored verbatim, -gpu num=n_gpu.

LSB_DJOB_NUMPROC then = gpu_cpu_cores (e.g. 8) -> mpirun -np 8 MMPBSA.py.MPI.
Cap still available via mmgbsa_n_cpu.

## Rendered headers (verified live)
rusage default  -> -n 8 + rusage[...,ngpus_physical=1] + -gpu num=1:j_exclusive=no
slots           -> -n 1 + -gpu num=1
n_slots=1 (FEP) -> -n 1 + -gpu num=1
gpu_cpu_cores=16-> -n 16 + ngpus_physical=1

## !!! CLUSTER-SPECIFIC VERIFICATION REQUIRED (cannot test from here) !!!
The EXACT GPU resource syntax LSF accepts is site-specific. `ngpus_physical`
is common, but your cluster may instead use:
  * `-gpu "num=1"` alone with `-n` already meaning cores (then use mode="slots"
    after confirming slots==cores), OR
  * a different resource name (e.g. `rusage[ngpus=1]` or `gpu` boolean).
ACTION: submit ONE job and check `bjobs -l <id>` / `bhosts -l` shows
  1*the GPU queue (NOT 8*) AND 8 cores granted. If bsub rejects the rusage line,
  switch gpu_request_mode to the syntax your admin confirms. Ask HPC support:
  "On the gpu queue, how do I request 1 GPU + 8 CPU cores in one job?"

## Tuning
  * gpu_cpu_cores: cores for scoring (8 default; 16 if node fair-share allows).
  * mmgbsa_n_cpu: cap MPI ranks below cores (memory pressure control).
  * Set gpu_cpu_cores=1 to revert to the final53 single-core behavior.

## Carried forward
final46/47/49/50/51/52/53.
