# v2.5.0 build final24: resume guard + batch delegation

## 1. Resume guard (run_pipeline)
Build + 5 ns MD is ~15 min/ligand on a GPU and is the expensive part. If the
workdir already has solvated.prmtop + solvated.inpcrd + production.nc, the
pipeline now SKIPS build+MD and scores only.
* Default: ON (resume). Re-runs after an analysis failure are near-instant.
* `--no-resume` forces a clean rebuild + fresh MD.
* Partial artifacts are handled: topo-only rebuilds MD; traj-only rebuilds topo
  then reuses the trajectory.

Verified (stubbed stages):
  clean+resume      -> build, md, score
  full+resume       -> score only          <-- the throughput win
  full+--no-resume  -> build, md, score
  topo-only+resume  -> build, md, score
  traj-only+resume  -> build, score (md reused)

## 2. run_batch delegates to run_pipeline
Per-ligand work now calls run_pipeline (clone of args -> lig_dir) instead of
duplicating build/MD/score. Benefits: the resume guard and EVERY analysis fix
(interval/comment, numpy.compat module, GB radii, IFBOX=0, MMPBSA diagnostics)
apply uniformly to batch runs. A failing ligand is logged and skipped; the batch
continues, then batch_aggregate ranks the survivors.

Verified (stubbed stages, 3 ligands incl. 1 designed to fail):
  round 1: build+md+score for all 3 (BAD raises, batch continues)
  round 2 (resume): score only for all 3
  batch_aggregate: parsed 3/3, ranked best->worst correctly, wrote
    binding_energies.tsv, binding_energies_ranked.md, INDEX.html, 3 reports.
