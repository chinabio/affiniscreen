# v2.5.3 build final49: make the MM-GBSA analysis kit / report actually fire

Context: the user's FINAL_RESULTS_MMGBSA.dat came from a MANUAL rescue command,
so its name (MMGBSA) + location differ from what the GUI pipeline writes
(<wd>/lig_<name>/mmgbsa/FINAL_RESULTS_MMPBSA.dat). Three latent bugs meant the
report/kit could silently no-op.

## Fix 1 -- analysis_kit shell scripts are now EXECUTABLE
run_analysis.sh / run_screen_analysis.sh were stored mode 600. submit.py's
analysis tail gates on `[ -x .../run_analysis.sh ]`, so the cpptraj
RMSD/RMSF/hbond/contacts kit was silently SKIPPED on Amber-engine runs. All
*.sh in the package are now packaged 0755 (verified in-zip).

## Fix 2 -- filename tolerance: FINAL_RESULTS_MM[PG]BSA.dat
mmpbsa_report.resolve_targets() and batch_aggregate now accept either
FINAL_RESULTS_MMPBSA.dat (pipeline) OR FINAL_RESULTS_MMGBSA.dat (older/manual),
including a bare .dat sitting in the workdir root. MMPBSA is preferred if both
exist. This makes the user's existing rescue file reportable AND future-proofs
against the naming variant.

## Fix 3 -- standalone report from any .dat (no MD rerun)
mmpbsa_report already exposes main()/__main__. Confirmed runnable as:
    python -m amber_md.mmpbsa_report /path/to/FINAL_RESULTS_MMGBSA.dat
    python -m amber_md.mmpbsa_report /path/to/workdir         # auto-find
-> writes FINAL_RESULTS.report.html next to the .dat. The is_file branch now
accepts ANY .dat name (covers _MMGBSA).

## NOT changed
* GUI command building, MD, mmgbsa.py output name (pipeline keeps writing the
  canonical FINAL_RESULTS_MMPBSA.dat under lig_<name>/mmgbsa/). The GUI path was
  already correct and self-contained (MMGBSAAnalyzer.run() generates the report
  inline); these fixes harden the edges + the manual/older-file case.

## Carried forward
final46 driver hardening (persisted resname + 0-atom mask abort),
v2.5.3 submit.py activate/abspath fix, final47 anti-storm guard.

## Verified
* changed .py compile; filename tolerance unit-checked; .sh exec bits confirmed
  0755 in the packaged zip.
