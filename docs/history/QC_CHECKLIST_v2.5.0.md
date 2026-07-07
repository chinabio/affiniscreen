# Amber MD / OpenFE Workflow v2.5.0 -- Click-through QC Checklist

Run top-to-bottom on the cluster (login node + `gpu` queue reachable).
Tick each box; note the build is v2.5.0 with legacy pages removed.

## 0. Launch
- [ ] `streamlit run amber_md/gui/Home.py` starts with no traceback.
- [ ] Landing page shows the v2.5.0 welcome table (5 pages only).
- [ ] Sidebar lists exactly: Setup & Launch, Job Monitor, Results - Single,
      Results - Compare, FEP Campaign (5 pages).
      No "Single Ligand", "Batch Screen", "ABFE Amber", "RBFE Amber",
      "ABFE/RBFE OpenFE", or old "Results" pages.

## 1. Setup & Launch -- MM-GBSA (Amber)
- [ ] Pick protein + a multi-record SDF; ligand records table populates.
- [ ] Method = MM-GBSA, Engine = Amber.
- [ ] Expand "MM-GBSA settings": set igb=5, salt=0.2, enable decomposition,
      set a decomp mask, set start/end frame + stride.
- [ ] Expand "Protonation (Amber)": untick auto, add `D:198:GLH` override.
- [ ] HPC expander shows **no** "Max concurrent GPU jobs" field; caption says
      jobs are submitted at once (LSF queues).
- [ ] Click **Validate** -> no errors. Click **Launch**.
- [ ] Generated command contains: `--config <path>`, `--salt 0.2`, `--decomp`,
      `--no-protonation`, `--protonate D:198:GLH`.
- [ ] Open the written `wizard_config.json` -> igb/decomp/frames/stride present.

## 2. Setup & Launch -- ABFE (Amber)
- [ ] Method = ABFE, Engine = Amber.
- [ ] Expand "Advanced Amber FEP": set a custom lambda schedule
      (must start 0.0, end 1.0), masks, temperature, cutoff.
- [ ] Boresch: leave "Auto-pick" ticked.
- [ ] Launch -> command routes through `python -m amber_md.fep_driver --mode abfe`
      with `--config` and `--auto-boresch`.
- [ ] Reject test: enter lambdas NOT starting at 0.0 -> validation error.

## 3. Setup & Launch -- RBFE (OpenMM/OpenFE)
- [ ] Method = RBFE, Engine = OpenMM; provide >=2 ligands.
- [ ] Single-ligand RBFE is blocked with a clear warning.
- [ ] Launch -> OpenFE planner command emitted; `openfe_network.yaml` written.

## 4. Job Monitor
- [ ] Active LSF jobs table renders (or "No active jobs").
- [ ] Tail-log expander works for a running job.

## 5. FEP Campaign -- Atom mapping (NEW)
- [ ] Point at a planned network dir.
- [ ] **Atom mapping** tab is first. With no artifacts -> info message.
- [ ] After planning (edges.csv present): edge table renders with
      core_heavy / perturbed_atoms / score columns.
- [ ] Perturbation-network Mermaid graph renders.
- [ ] "Per-edge MCS SMARTS & masks" expander shows SMARTS per edge.
- [ ] "Run edges" tab: submitting shows NO "GPU slots full" wait message;
      all edge x repeat jobs submit immediately.

## 6. Results - Single Molecule
- [ ] Auto-targets last launched ligand (or pick a dir).
- [ ] Headline DG_bind + per-leg breakdown + convergence render.
- [ ] MM-GBSA per-residue decomposition shows when present.

## 7. Results - Compare & Rank + Promote-to-FEP (NEW)
- [ ] Point at a batch parent dir; ranking table + chart render.
- [ ] Paste experimental DG -> RMSE / scatter vs experiment appears.
- [ ] **Promote MM-GBSA hits -> FEP** section lists hits lacking a FEP scaffold.
- [ ] Choose top-N, mode=relative, set masks; submit.
- [ ] Each `fep_<ligand>/` dir is created with `submit_fep.sh`;
      status table shows OK / SKIP / FAIL correctly.
- [ ] mode=absolute -> command uses `--mode abfe --auto-boresch`.

## 8. Regression guards
- [ ] Grep the tree: no `MAX_GPU = 8`, no `while count_my_jobs`,
      no `while slots_fn`, no "Max concurrent GPU jobs".
- [ ] `WorkflowConfig.save()`/`load()` round-trips `auto_protonation`
      and `protonation_overrides`.
- [ ] Every `.py` carries exactly one author header (no duplicates).

---
_Generated for v2.5.0. Legacy pages 2,3,4,5,6,7,8 removed; their functionality
lives in Setup & Launch, FEP Campaign, and the two Results pages._
