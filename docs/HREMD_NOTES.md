# HREMD in amber_md — architecture notes and why it stays opt-in

**Status:** `--hremd` is implemented but **deliberately opt-in and NOT
recommended on the the GPU queue queue.** This document records why, so the decision
is not silently reversed later.

## TL;DR
The Amber HREMD path in this package (`build_lsf_hremd`, fep.py) uses Amber's
native **synchronous** multisander replica exchange:

```bash
mpirun -np {nwin} pmemd.cuda.MPI -ng {nwin} \
       -groupfile groupfile.prod -rem 3 -remlog rem.log
```

`-np N` launches **N MPI ranks**, `-ng N` declares **N replica groups** —
**one persistent MPI rank (and one live CUDA context) per replica.** All ranks
step **in lockstep** and exchange via an MPI collective every `numexchg` block.

**Consequence: the number of GPUs required equals the number of windows.**
A 44-window `complex_vdw` leg needs ~44 GPUs **concurrently on one host**
(`#BSUB -R "span[hosts=1]"`). On the GPU queue (~8 GPUs/host) this cannot run.

## "But OpenFE runs replica exchange on 8 GPUs — why can't Amber?"
Because OpenFE/openmmtools use a fundamentally different, **asynchronous**
replica-exchange architecture:

| | OpenFE / openmmtools | Amber `pmemd.cuda.MPI -ng` (this package) |
|---|---|---|
| REX model | Asynchronous, Python-driven | Synchronous, MPI lockstep |
| Engine | OpenMM, in-process swappable contexts | pmemd.cuda.MPI, persistent per-rank |
| Exchange | Python/numba math on energies (`_mix_replicas`) | MPI collective between live ranks |
| GPUs needed | Decoupled — 44 states on 1–8 GPUs (time-sliced) | = #replicas (44 wants ~44 GPUs) |
| 44 windows on 8-GPU queue | Yes (round-robin onto available GPUs) | No (needs 44 concurrent, or oversubscribe → OOM / no gain) |

`openmmtools.multistate.ReplicaExchangeSampler` is **one process** that holds all
thermodynamic states, propagates them one/few at a time onto whatever GPUs exist,
then mixes states with pure CPU math. GPU count is a throughput knob completely
decoupled from replica count. Amber pmemd has **no equivalent async scheduler** —
its only multi-GPU mode is the lockstep groupfile model above. Oversubscribing 44
ranks onto 8 GPUs via MPS runs all 44 contexts concurrently (each a full ~113k-atom
system, ~1.5 GB+), competing for memory/compute — expect OOM and no speedup.

**Therefore Amber cannot replicate OpenFE-style "REX on a small GPU pool." It is
an engine-level limitation, not a config flag or a bug in this package.** To get
async few-GPU replica exchange, use OpenFE/openmmtools (or GROMACS `-multidir`
HREX), not Amber.

## Known issues in the current `build_lsf_hremd` (if it is ever enabled)
1. **Window-count bug:** uses `len(self.cfg.lambdas)` (the generic ~11-window
   set) instead of `self._active_lambdas` (the per-stage 21/44 schedule used by
   `build_lsf_array`). For two-stage ABFE this builds the wrong number of
   replicas. Must be reconciled before use.
2. **Not wired into the GUI** — `--hremd` / `--exchange-freq` are CLI-only; the
   GUI launch builder does not emit them. Ticking anything in the GUI will not
   enable HREMD.
3. **Resume / analyze paths** were written for the independent-array model and
   have not been validated against a coupled REX job.

## Recommendation
- **Use the plain-MD independent array** (`build_lsf_array`, the default). It is
  the correct architecture for an 8-GPU queue: embarrassingly parallel, degrades
  gracefully (N run, rest pend), and is per-window debuggable.
- **Do not enable `--hremd`** unless you can allocate `nwin` GPUs on a single
  host (not the case on the GPU queue for the 44-window vdw leg) **and** items 1–3
  above have been fixed and tested.
- If replica-exchange sampling quality matters for a specific target, run that
  target in **OpenFE** (async REX) rather than forcing Amber into a mode its
  engine does not support.

## References
- Amber: `pmemd.cuda.MPI` multi-GPU is only for "methods requiring multiple
  simulations to communicate" (TI, REMD); single-sim multi-GPU scaling is poor.
- openmmtools `ReplicaExchangeSampler`: one-process async multistate sampler;
  `_mix_all_replicas_numba` performs swaps as CPU energy math.
