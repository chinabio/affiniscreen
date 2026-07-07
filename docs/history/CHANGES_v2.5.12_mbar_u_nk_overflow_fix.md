# v2.5.12 (final68) — MBAR/BAR u_nk soft-core overflow fix

## Symptom
The first full ABFE cycle ran end-to-end but produced an UNTRUSTED result with
`dG_bind = +55.4 kcal/mol`. Per-leg `dG_estimators.csv` showed TI sane but
MBAR/BAR wildly wrong on EVERY leg — including the fully-complete ones:

| leg                | TI      | BAR     | MBAR      |
|--------------------|---------|---------|-----------|
| complex_decharge   | -30.59  | +58.17  | **-1148.36** |
| solvent_decharge   |  +8.89  | +33.87  |  **-114.77** |
| solvent_vdw (23/23)|  -2.46  |  -6.91  |  **-334.26** |

A decharge leg cannot be -1148 kcal/mol. Because TI was correct, this was an
ANALYSIS bug, not a sampling/physics problem — no re-run could fix it.

## Root cause (confirmed on real prod.out data)
The alchemlyb AMBER parser builds the reduced-potential matrix `u_nk` from the
`MBAR Energy analysis` blocks. At a soft-core window, re-evaluating a coupled
config at a (near-)decoupled end state produces an atom-overlap singularity:
individual `u_nk` CELLS blow up to ~1e6 kT (observed: `u_nk[lambda=0.0]=2.4e6`
for a lambda=0.5 sample; 8 such frames out of 2500 at lambda=0.5 alone).

The final60 guard `_sanitize_u_nk` tried to CLAMP each cell to
`per-row median +/- _UNK_SANE_MAX_KT`. But when several cells in a frame
overflow, the per-row median is itself enormous, so the cap never bites and the
1e6 values flow into pymbar. pymbar's solver then fails
(`DLASCL parameter 4 had an illegal value`, `min nearest-neighbour overlap
0.0000`) and returns free energies inflated by ~the number of states. TI was
spared because it only uses `dHdl` (sanitized separately).

The pre-final60 approach (drop every row with any bad cell) was also wrong: at a
fully-overflowing high-lambda window it emptied an entire lambda group, making
the matrix rectangular -> `Shape of passed values is (n,n), indices imply
(n+1,n+1)`.

## Fix
Rewrote `FEPAnalyzer._sanitize_u_nk` to DROP whole frames whose reduced
potential is non-finite or exceeds `_UNK_SANE_MAX_KT` (1e4) kT at ANY state,
done PER lambda group so no window is emptied. If dropping would leave a group
with fewer than `_UNK_MIN_KEEP` (50) frames, the group instead keeps its
least-extreme 50 frames, so the reduced matrix stays square and every state
stays populated. TI is unaffected.

New constant: `_UNK_MIN_KEEP = 50`.

## Verification (on the user's real solvent_vdw prod.out)
Through the actual patched method:

| lambda | raw max u_nk | dropped | survivors | clean max | cols | own-state mean |
|--------|--------------|---------|-----------|-----------|------|----------------|
| 0.000  | 0.0e0        | 0       | 2500      | 0.0 kT    | 23   | 0.0 |
| 0.500  | 2.43e6       | 26      | 2474      | 9948 kT   | 23   | 0.0 |

The cleaned lambda=0.5 profile is a smooth convex bowl with its minimum at the
own state — well-formed MBAR input. Matrix shape preserved (23 columns).

NOTE: a true 23-state MBAR free energy could not be recomputed in-sandbox (only
2 of 23 windows were available); validation confirms the u_nk matrix is now
physically well-formed and overflow-free, which is the precondition pymbar
needs. The headline number should be re-confirmed once complex_vdw completes
under 2.5.11's schedule and the full leg is re-analyzed.

## Release gate
- 142 entries, no files dropped, no unintended changes.
- byte-compile: ALL OK; version 2.5.12 (final68).
- bash -n (array + HREMD generated scripts): clean.
- pytest: 29 passed.

## Lineage
  * exit 141 (SIGPIPE eq gate)              -> v2.5.7
  * host / outage kills                      -> v2.5.8
  * exit 2 (duplicated if -> syntax)         -> v2.5.9
  * exit 255 (box-drift grep wrong file)     -> v2.5.10
  * exit 71 / residual 255 (lambda spacing)  -> v2.5.11
  * MBAR/BAR garbage (u_nk cell overflow)    -> v2.5.12 (this)
