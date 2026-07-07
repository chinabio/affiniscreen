CHANGES_abfe_submit_flag (v2.5.3 / final56)
===========================================
Date: 2026-06-09

SYMPTOM
-------
After final55 fixed the --config crash, the Amber ABFE driver ran to completion
(topology, Boresch, 76 windows, cycle-closer script) but NO jobs appeared in the
LSF queue. launch_abfe_0.log showed "wrote N lambda windows ..." for every leg
and ended at "Wrote cycle-closer script: ...", with no "submitted MD=..." lines.

ROOT CAUSE
----------
fep_driver.run_fep() only bsub's windows when args.submit (or args.resume) is
truthy; otherwise it takes the else-branch and merely writes the per-window LSF
scripts (build_lsf_array / _build_analyze_lsf) without calling bsub. --submit is
a bare store_true flag on the driver.

The GUI ABFE/Amber command builder never appended --submit, even when the user
selected the "GPU queue (LSF)" execution target. (The MM-GBSA/OpenMM branch
already gated `--submit gpu` on that same toggle; ABFE was simply never wired.)

FIX
---
0_Setup_and_Launch.py:
  - ABFE/Amber branch: append "--submit" when
    exec_target startswith "GPU queue" (mirrors the MM-GBSA gate). On
    "Local host" it is omitted, preserving dry-prep/scripts-only behaviour.
  - Launch UI: _is_queue now also recognises ABFE/Amber + GPU queue, so the
    user sees the "Submitting to LSF GPU queue" confirmation (and job-id
    capture) instead of a bare "Launched (pid ...)".

IMPACT
------
- Amber ABFE on "GPU queue (LSF)" now actually submits all four legs and the
  cycle-closer; ABFE_RESULT.txt will be produced when both legs finish.
- "Local host" behaviour unchanged (prep-only).
- No change to MM-GBSA / RBFE / OpenFE paths.

FILES CHANGED
-------------
- amber_md/gui/pages/0_Setup_and_Launch.py
