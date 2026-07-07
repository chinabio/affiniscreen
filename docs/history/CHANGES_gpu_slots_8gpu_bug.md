# v2.5.3 build final53: FIX — Amber MM-GBSA reserved 8 GPUs (8*the GPU queue)

## Symptom (user bjobs)
  4688644 ... RUN gpu ... 8*gpu-node-02 amberMD
A single-GPU pmemd.cuda MM-GBSA job was holding 8 GPUs: 7 idle, slow to
schedule, and it did NOT speed up the run (still ~13 ns/day).

## Root cause
submit.py._header() GPU branch set `#BSUB -n` from n_cpu (=8), intending
"8 CPU cores + 1 GPU on one host" (with `-gpu num=1`). But on the the cluster `gpu`
queue, LSF SLOTS map to GPUs, so `-n 8` reserved 8 GPUs regardless of
`-gpu num=1`. The 8-core intent (for the in-job MM/GBSA MPI step) was wrong for
this slot-to-GPU mapping.

## Fix
* config.py HPCConfig: new field `n_gpu_slots: int = 1`.
* submit.py GPU branch: `#BSUB -n` now uses n_gpu_slots (=1), NOT n_cpu. An
  explicit n_slots arg still overrides (FEP windows unaffected).
* The MM/GBSA MPI step reads $LSB_DJOB_NUMPROC at runtime, so it simply uses
  whatever LSF grants (now 1) — no breakage; MM/GBSA scoring is fast anyway.

## Effect
GPU header now: `#BSUB -n 1` + `#BSUB -gpu "num=1"` -> 1*the GPU queue.
Frees 7 GPUs per job, schedules faster, identical run speed.

## Note on speed
This does NOT change ns/day (Amber MD was always using only 1 GPU). It stops
the WASTE of 7 reserved-idle GPUs and the scheduling penalty. The ~13 ns/day
seen earlier was likely GPU contention with the concurrent 50 ns job; read the
true rate from `grep ns/day jobs/prod.out` on a clean run.

## If you WANT more cores for MM/GBSA scoring
Set n_gpu_slots higher ONLY if your site maps GPU-queue slots to cores, not
GPUs. On the cluster `gpu`, keep it at 1.

## Verified live (imported package, rendered header)
* default GPU header -> `#BSUB -n 1`, `-gpu num=1`.
* explicit n_slots=2 still emits -n 2 (FEP override intact).

## Carried forward
final46 driver hardening, v2.5.3 submit activate fix, final47 anti-storm guard,
final49 report-kit, final50 unified report button, final51 prod-ns unification,
final52 prod-ns config plumbing (GUI 10 ns now honored).
