# v2.5.3 build final55: FIX bsub rejection -- GPU specified twice

## Symptom (live bsub, exit 255)
  STDERR: GPU number in GPU resource requirements and -gpu cannot be specified
          at the same time. Job not submitted.

## Cause
final54 rusage mode emitted BOTH rusage[ngpus_physical=1] AND -gpu "num=1".
This LSF allows only ONE GPU specification.

## Fix (submit.py rusage branch)
Removed the #BSUB -gpu line in rusage mode. GPU requested ONLY via
rusage[ngpus_physical={n_gpu}]; CPU cores via -n {gpu_cpu_cores}.

## All three paths now emit EXACTLY ONE GPU spec (verified live)
  * rusage default       : -n 8 + rusage[mem,ngpus_physical=1]  (no -gpu)
  * explicit n_slots(FEP): -n N + -gpu "num=1"
  * slots legacy         : -n 1 + -gpu "num=1"

## Expected
  bjobs -w  -> 1*the GPU queue with 8 cores
  mmgbsa.log -> [MM/GBSA] using mpirun -np 8 MMPBSA.py.MPI (LSB_DJOB_NUMPROC=8)

## If still bounces
ask HPC exact resource name; fallback rusage[ngpus=1] or gpu_request_mode=slots.

## Carried forward
final46/47/49/50/51/52/53/54.
