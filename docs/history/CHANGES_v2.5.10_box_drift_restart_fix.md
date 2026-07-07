# v2.5.10 (final66) — FIX: prod box-drift auto-restart never triggered

## Symptom
A clean v2.5.9 `complex_decharge` run got **13/15 windows to full completion**
(eq -> prod -> "all stages completed OK") — the eq gate is fixed. But 2 windows
(lambda 0.275, 0.725) failed with **exit 255** during `prod`, and the log said:

```
[run_stage] lambda=0.275 prod attempt 1/5 from eq.rst
ERROR: Calculation halted.  Periodic box dimensions have changed too much ...
STAGE FAILED: ... stage=prod rc=255 (non-box error)
```

i.e. a textbook box-drift halt was misclassified as a "non-box error" and the
window died on attempt 1 — the auto-restart loop never engaged.

## Root cause
pmemd.cuda writes the box-drift halt message
(`Periodic box dimensions have changed too much from their initial values`) to
**STDOUT**, which LSF captures in the job `.out`. It is **not** written into
`prod.out` (the `-o` mdout, which just ends at the last energy frame).

`run_prod_with_restart()` detected drift with:

```bash
if grep -qiE "...changed too much..." prod.out 2>/dev/null; then  # never matches
```

Because `prod.out` never contains the string, the grep always failed, the retry
branch was skipped, and execution fell through to the generic "non-box error"
path -> `return 255`.

## Fix
Capture pmemd's stdout+stderr and grep that:

```bash
pmemd.cuda -O -i prod.in ... -o prod.out -r prod.rst -x prod.nc \
           > prod.console.${attempt} 2>&1
rc=$?
cat prod.console.${attempt} || true        # keep it visible in the LSF log
...
if grep -qiE 'box dimensions have changed too much|changed too much from their initial|Periodic box dimensions have changed' \
        prod.console.${attempt} prod.out 2>/dev/null; then
    # restart from latest rst (prod.rst if present, else eq.rst)
    attempt=$((attempt+1)); continue
fi
```

The fallback diagnostic grep now also scans `prod.console.<attempt>` so a
genuine non-box 255 (NaN/instability) still surfaces its cause.

This affects the **array builder only**. The HREMD builder runs production as a
single coupled `mpirun -rem 3` replica-exchange job and has no per-window prod
restart by design, so it is unaffected.

## Verification
- `bash -n` clean for both array + HREMD generated scripts (release gate).
- Functional simulation with a fake `pmemd.cuda` that emits the real Amber
  box-drift message to STDOUT on attempt 1 then succeeds: the loop restarts and
  reports `prod completed OK (attempt 2)` — exit 0. (Old code would have exited
  255 on attempt 1.)
- Full unit suite: 29 passed.

## Note on the 2 failed windows from the v2.5.9 run
Those windows are physically recoverable — box drift is a known NPT/GPU-grid
issue that a restart fixes. Re-run them (the self-heal `--resume` will relaunch
only the unfinished windows), and with this fix the in-job restart loop will
ride through the drift automatically.

## Lineage (all distinct failure modes, now fixed)
  * exit 141 (SIGPIPE, eq gate `head` under pipefail)     -> v2.5.7
  * host-specific / outage kills                           -> v2.5.8
  * exit 2 (duplicated `if` -> bash syntax error)          -> v2.5.9
  * exit 255 (box-drift restart grepped wrong file)        -> v2.5.10 (this)
