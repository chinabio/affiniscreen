#   2.5.77 - FINAL clean-run build. Fixes from the 5 ns diagnostic run:
#            (B) mbar_lambda %.3f truncated 4-dp lambdas -> 11/30 windows no
#            u_nk -> MBAR dead: schedules now 3-dp + writer %.4f. (A) GUI prod
#            slider default 5 ns overrode 20 ns config -> GUI defaults ->20 ns.
#            (C) data-driven inserts 0.656/0.706 at the decharge dV/dl well
#            (dV/dl +49->-93, overflows at ~0.688); decharge 30->32 windows.
#            (D) preflight_abfe.check_prodin_consistency() guards both bugs.
#   2.5.75 - CONVERGENCE-GUARANTEED (10-14 day budget): nstlim_prod 10->20 ns.
#            160 win x 20 ns = 3200 ns ~7.8 d wall @8 GPU. Schedules unchanged
#            (decharge 30 / dense vdw 50 already adequately spaced). Doubling
#            ns targets the fwd/rev ~16 kcal/mol drift directly.
#   2.5.74 - SAMPLING BUDGET (8 GPU x 1 week): nstlim_prod 5->10 ns;
#            decharge 26->30 win (well <=0.0375); dense vdw 40->50 win
#            (danger zone <=0.0125). 146 win x 10 ns = 1460 ns ~3.8 d wall.
#   2.5.73 - CONVERGENCE FIXES (all legs): nstlim_prod 2->5 ns; decharge
#            lambdas densified 0.725-0.925 (21->26 win, MBAR overlap); vdw
#            dense 40-win default ON; convergence_analysis.py None-guard
#            (lambda=1.0 crash); analyze LSF now RUNS convergence_analysis.py
#            (writes convergence.csv/.png); HREMD _active_lambdas fix.
#   2.5.72 - DOCS: docs/HREMD_NOTES.md -- records that Amber pmemd -ng HREMD
#            is SYNCHRONOUS (needs ~1 GPU/replica on one host), so it cannot
#            do OpenFE-style few-GPU async REX; HREMD stays opt-in. No
#            behaviour change.
#   2.5.71 - PERF (box buffer 16->12 A; OpenFE-aligned, root-cause for 16 A
#            fixed in 2.5.70) + ships tools/convergence_analysis.py. HREMD
#            kept opt-in. See CHANGELOG_v2.5.71.md.
#   2.5.70 - FIX (complex_decharge protocol): hard-vdW decharge (stage-aware
#            gti_vdw_sc=0 / gti_chg_keep=0) kills the charge-on-softcore
#            singularity; decharge schedule densified 15->21 windows; MD-only
#            NaN/**** prod health guard. See CHANGELOG_v2.5.70.md.
#   2.5.61 - RELEASE (timestep stability + restraint-leg analysis closure).
#            FUNCTIONAL FIX (all FEP legs): production dt 0.002 -> 0.001 ps.
#            dt=0.002 proven (multi-seed A/B) on the pmemd.cuda stability edge
#            for lig_12944901: restraint mid-band (lambda 0.10-0.40) crashed 3/3
#            seeds; complex_decharge lambda_0.175 also died. dt=0.001 removes the
#            crash class (2x wall-clock per ns).
#            GUI FIX: ABFE launcher never passed --dt (dt from FEPConfig.dt_ps)
#            yet computed --nstlim-prod with a hard-coded /2 (2 fs) -> decoupled;
#            changing dt alone would HALVE simulated time. Fixed BOTH.
#            FINDING 1: FEPAnalyzer cannot analyze the Option-A restraint leg
#            (needs icfe=1/clambda; leg is icfe=0) -> 0/N windows.
#            FINDING 2: Option-A leg is not a sound numerical-TI source
#            (1/lambda endpoint singularity). Restraint = ANALYTIC -11.45.
#            DOCS: docs/release_v2.5.61/ bundles the crash-investigation handoff
#            and the apply_dt001_* patch scripts. VERSION TRACK resynced (VERSION
#            + README were 2.5.51a). See CHANGES_v2.5.61_timestep_and_restraint_analysis.md.
#   2.5.60 - FIX (restraint-leg recovery ladder): mid-band windows (lambda
#            0.10-0.32) hit pmemd.cuda "box changed too much" repeatedly and
#            EXHAUSTED the regrid ladder (10 restarts -> rc=255) despite a
#            perfectly stable box (density ~1.016, no clash). Root cause: once
#            _drift_count>=3 the ladder FROZE in tier-2 and restarted from
#            orig_eq.rst (pristine origin) every time, regenerating the same
#            origin grid that tripped the check -> zero forward progress.
#            Amber docs: the error fires even on a healthy box; remedy is to
#            restart FROM THE LATEST restart file. FIX: classify each drift as
#            CLASH (density outside [0.85,1.15], |PRESS|>5000, or NaN VDWAALS)
#            vs BENIGN. Benign -> always regrid from $good_rst (advance). Clash
#            -> rewind to orig_eq.rst AT MOST ONCE (_clash_rewound), then
#            resume from $good_rst. Preserves the lambda=0.150 clash handling
#            while unblocking benign-fluctuation windows.
#   2.5.59 - FIX (restraint-leg prod instability): the restraint production
#            mdin was the ONLY prod stage WITHOUT a vlimit guard. 4 of 14
#            windows in the lambda 0.10-0.26 band detonated mid-run (single
#            2 ps block: BOND->1e8, VOLUME -20%, TEMP=NaN; energies healthy
#            in the prior block -> instantaneous clash, not drift/timestep).
#            FIX1: _prod_in_restraint now emits vlimit (cfg.vlimit floored at
#            20.0) for parity with every other prod stage. dt and lambda
#            schedule unchanged (diagnostic showed neither was the cause).
#            FIX2: run_prod_with_restart now NaN-gates a rc=0 prod.out --
#            a detonated-but-exit-0 window is failed (rc=255) instead of
#            silently poisoning MBAR/TI.
#   2.5.58 - CRITICAL FIX (regression in 2.5.54): skinnb was emitted inside
#            the &cntrl namelist, which pmemd does NOT accept -> every window
#            aborted at the dens stage with "Cannot match namelist object name
#            skinnb" (rc=2), producing 0/N windows and no dG. skinnb is an
#            &ewald variable (verified: AMBER devs list + ParmEd ewald namelist,
#            pmemd default 2.0). FIX: skinnb=3.0 now emitted as a separate
#            &ewald block after &cntrl in eq/dens/prod/restraint-prod mdins AND
#            the recovery mdin. All four mdin renders validated (cntrl+ewald
#            separate, both namelists closed, skinnb in &ewald only). _min_in
#            left unchanged (no PME dynamics). ANYONE ON 2.5.54-2.5.57 MUST
#            UPGRADE: those versions cannot run a single MD step.
#   2.5.57 - DOCS: PACKAGE_COMPARISON.md -> v1.1. BAT.py now GROUNDED (real 7MB
#            source uploaded). Read from code: routes express/dd/sdr/dd-rest
#            (fe_type), components m/n/c/e/v/r, attach_rest=10 weight windows,
#            decoupling lambdas=23 (express-MBAR), TI-Gaussian or MBAR, MC
#            barostat (barostat=2), HMR dt=0.004, cut=9, 298.15K, multi-DOF
#            restraints (dist/angle/dih/COM). All 4 references now permanently
#            captured -> no re-upload needed. No code change.
#   2.5.56 - DOCS: added PACKAGE_COMPARISON.md - grounded side-by-side of ALL
#            reference ABFE codebases so they need not be re-uploaded. Read from
#            source: FEP-SPell-ABFE (taup=2.0, scalpha=0.2/scbeta=50, cut=9,
#            5ns legs, vdw=44/restraint=16 windows -- this workflow's lineage);
#            the reference platform (OpenMM, MC barostat npt_mc_freq=25, 3-vector lambda, recipes
#            v44/e29/r01, Boresch k=25); FEP+ (108/60 windows, 10ns, REST+GCMC).
#            BAT.py NOT analyzed: uploaded zip is a 96-byte stub (flagged in doc).
#            No code change.
#   2.5.55 - DOCS: added COMPARATIVE_ANALYSIS.md (design rationale + grounded
#            cross-tool comparison). FEP+ protocol extracted from the project
#            .msj files: complex leg = 108 lambda windows, solvent leg = 60,
#            production = 10 ns/window (time=10000ps), REST + lambda-hopping
#            (HREX) + GCMC water (muVT), scale_solvent_vdw=0.75; full multisim
#            wall time 12h24m / ligand. Documents WHY the pmemd "box changed
#            too much" failure (fixed nonbond grid) cannot occur in Desmond/
#            OpenMM (FEP+/the reference platform), confirming the v2.5.54 skinnb=3.0 fix as the
#            engine-specific root cause. No code change.
#   2.5.54 - ROOT-CAUSE FIX for the "box changed too much" false-fail. Full
#            21-window restraint-leg triage proved the box was NOT drifting:
#            8/21 windows died at 0.1-22% of prod with PRESS bounded, density
#            steady ~1.016, volume <0.2% variation, NO NaN; prod.out truncated
#            mid-run with NO error in it -- pmemd.cuda printed the halt ONLY to
#            the console. eq(500k)+dens(250k) end-states identical for passing
#            and failing windows. Cause: nonbond cell-list/skin error, not
#            drift -- skinnb unset (pmemd.cuda default 2.0 A) with cut=10.0,
#            too tight on the large complex box. FIX: skinnb=3.0 added to prod
#            (TI+restraint), eq, dens AND the regrid recovery mdin (which also
#            lacked it -- why regrids kept re-tripping). 2.5.52 taup=5/barostat=1
#            kept (complementary).
#   2.5.53 - Hardened the 2.5.52 box-drift "accept at completion" gate. Step
#            count (>=90% nstlim) is now necessary but NOT sufficient: the run
#            must ALSO pass a health check (no NaN in prod.out AND last Density
#            in [0.85,1.15] g/cm^3) before being accepted without rewind. A run
#            that reached 90% but is physically degraded continues the recovery
#            ladder instead of being falsely stamped converged.
#   2.5.52 - BOX-DRIFT FALSE-FAIL FIX (lambda=0.100 case). Two changes:
#            (1) prod.in (restraint + TI legs) now uses taup=5.0, barostat=1
#                (gentle Berendsen) instead of taup=2.0. Matches the already-
#                working dens.in/eq.in; stops pmemd.cuda fixed-grid "box changed
#                too much" halts caused by NPT volume excursions. the reference platform (OpenMM,
#                MC barostat) has no such failure mode because OpenMM auto-
#                regrids; on pmemd a gentle Berentsen coupling is the safe analog.
#            (2) run_prod_with_restart() no longer REWINDS a healthy run. If a
#                box-drift halt occurs at >=90% of nstlim, the segment is accepted
#                as converged (return 0) instead of restarting the full segment
#                from settled coords (which re-drifts to the same halt -> the
#                observed 10x loop + false STAGE FAILED). Mirrors the reference platform resume-
#                from-checkpoint design (never restarts a whole segment).
"""AffiniScreen."""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
# Portions ported/adapted from FEP-SPell-ABFE (freeenergylab, MIT License).

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

# ---------------------------------------------------------------------------
# SINGLE SOURCE OF TRUTH for the package version.
# Bump this on EVERY change/fix. It is shown on the GUI homepage and stamped
# into every generated LSF script header (see amber_md.version.lsf_banner).
#
# Changelog (newest first):
#   2.6.0  - GUI STREAMLINED. Removed Amber ABFE and Amber RBFE from the
#            Streamlit GUI (never validated end-to-end); the underlying
#            engine code (fep_driver, rbfe_map, abfe_* modules) is retained
#            for programmatic/CLI use. FEP Campaign is OpenFE-only; the
#            'Promote MM-GBSA -> Amber FEP' action was removed. Repository
#            reorganized for GitHub (LICENSE, docs/, examples/, history/).
#   2.5.33 - RESTRAINT LEG, FEP-SPell-FAITHFUL (approach A). Reverts the v2.5.32
#            single-topology force-constant scaling (which exploded the box at the
#            dens stage: Density 0.126, illegal memory access on all 16 windows).
#            Root cause: build_restraint_topology RE-SOLVATED the dual-copy system
#            from scratch -> un-equilibrated water packing the barostat crushed.
#            FIX (per FEP-SPell abfe_alchemy_morph): build the dual-copy restraint
#            topology WITHOUT re-solvation -- carry over the EQUILIBRATED waters +
#            box (strip only the ligand, combine {LIG LIG protein_solv}, setBox +
#            ChBox to restore the oct box). Restraint leg is dual-copy :1/:2 TI
#            with FIXED-k Boresch on :1, ifsc=0, NO crgmask/reion (all FEP-SPell
#            settings). Keeps the 2.5.31m box-drift recovery.
#   2.5.32 - RESTRAINT LEG REDESIGN. The dual-copy (:1 real / :2 dummy) restraint
#            topology was annihilating a fully-charged ligand (restraint_crgmask=''
#            left :2 interacting), giving dV/dl -> -1174 at lambda=1 and a bogus
#            complex_restraint dG = -125 kcal/mol that corrupted dG_bind (MBAR
#            rank-deficient, BAR diverged). FIX: retire the dual-copy topology;
#            the restraint leg now runs on the REAL complex prmtop with IDENTICAL
#            fully-interacting TI ends and scales the Boresch force constants
#            rk(lambda)=lambda*k (0->full). decharge/vdw keep FIXED full k
#            (restraint held on). This is the standard ABFE restraint-introduction
#            leg; its dG is small and pairs with the analytic standard-state term.
#   2.5.31m- ROBUST box-drift recovery for EVERY leg. 2.5.31l regrid still failed
#            rc=1 because it restarted from the CORRUPT drift-halted prod.rst with
#            irest=1/ntx=5. Now: regrid restarts from the LAST KNOWN-GOOD frame,
#            reads COORDS-ONLY (irest=0,ntx=1), and DERIVES prod_regrid.in from
#            each leg's own prod.in via sed so masks/icfe/crgmask/clambda match.
#            Retries 5->8 with a coords-only final fallback.
#   2.5.31l- Box-drift REGRID step was DEAD: prod_regrid.in heredoc used a QUOTED
#            delimiter (<<'REGRID') so the intended runtime shell var ${LAM} was
#            written literally -> "clambda=${LAM}" -> pmemd "namelist not
#            terminated" -> all 5 regrid attempts failed -> box-drift windows
#            (e.g. lambda=0.750) never produced a real prod.nc. FIX: unquote the
#            heredoc delimiter so the shell expands ${LAM}; added a guard (exit 74)
#            that aborts if clambda is still non-numeric.
#   2.5.31k- FAIL FAST: added a PRE-FLIGHT Boresch geometry gate right after min
#            (validates min.rst). A broken restraint is now caught in seconds
#            (exit 73) instead of ~50 min later at the post-eq gate. Same
#            validator; mirrors the post-eq gate's import-robust heredoc.
#   2.5.31j- Dual-copy restraint topology placed the LIGAND OUTSIDE THE POCKET.
#            build_restraint_topology loaded the ligand from the standalone
#            LIG.mol2 frame while the protein came from the complex frame, so
#            tleap `combine` left the ligand ~63 A away -> RESTRAINT Dist ~39,000
#            even after 31i fixed the indices (BOND/VDW were already normal). FIX:
#            extract the ligand BOUND POSE from the complex via cpptraj and
#            transplant those coords into LIG.mol2 (params/types kept) before
#            combine. Added a post-build pocket sanity guard.
#   2.5.31i- TRUE ROOT CAUSE of the 49,000 kcal/mol restraint (it survived the
#            31e anchor-swap fix): the dual-copy Boresch index remap shifted
#            PROTEIN atoms by +n_lig instead of +2*n_lig. The dual-copy layout is
#            :1 lig / :2 lig / protein, and the original complex is protein-first
#            ligand-last, so two ligand copies precede the protein -> shift must
#            be 2*n_lig. The old +n_lig left every receptor anchor 37 atoms low
#            (wrong residue ~65 A away) -> boresch.RST iat=5,10567 (should be
#            10604) -> RESTRAINT ~49,000 -> gate FAIL r=65.21. Fixed shift +
#            added a guard rejecting any anchor that lands in the :2 dummy region.
#   2.5.31h- NumPy 2.0 vs Amber-22 parmed: pip install -e . (v2.5.31g) pulled
#            numpy 2.x into the amber-md env, which made Amber's parmed crash in
#            the post-eq Boresch gate (np.array(copy=False) -> ValueError),
#            rc=1 -> production refused. FIX: pin numpy<2 in pyproject AND make
#            validate_masks tolerate the error. min/heat/dens/eq all ran clean on
#            GPU; this dependency clash was the only blocker.
#   2.5.31g- PACKAGING: add pyproject.toml so `pip install -e .` works (was failing
#            with "neither setup.py nor pyproject.toml found"). Dynamic version read
#            from amber_md.__version__ (never drifts). parmed/pytraj deliberately
#            excluded from deps (Amber provides them via PYTHONPATH); MDAnalysis/
#            mdtraj/pymbar/alchemlyb under [analysis] extra, streamlit/nglview under
#            [gui]. Console scripts: amber-fep-driver, run-amber.
#   2.5.31f- RESTORE mdin_validator.format_issues(). It was orphaned as dead code
#            after the return in validate_paths(), so it was never defined; fep_driver
#            and tools/validate_mdin.py could not import it and the mdin content/-ref
#            gate was SILENTLY SKIPPED ("could not run ...; continuing") -- meaning
#            the TEMP0/&wt safety net we added did not actually run at submit time.
#            Gate is functional again. (v2.5.31e Boresch anchor-swap fix confirmed
#            working in production: correction now -11.4 kcal/mol, geometry PASS.)
#   2.5.31e- TRUE ROOT CAUSE of the 50,000 kcal/mol Boresch restraint: 
#            abfe_integration.boresch_dict_from_prep() SWAPPED ligand and receptor
#            anchors -- it wrote ligand indices into aA/bA/cA (receptor slots) and
#            receptor indices into A/B/C (ligand slots). _write_boresch_RST then
#            restrained the WRONG pair (~65 A apart) to the correctly-measured
#            r0 (~5 A) -> RESTRAINT ~50,000 kcal/mol -> GPU prod step-1 abort.
#            Fixed the mapping to canonical A/B/C=ligand, aA/bA/cA=receptor. The
#            v2.5.31d write-time guard remains as defense-in-depth.
#   2.5.31d- PROD STAGE FIXES (evidence from abfe_20260617_150311 lig_12944901):
#            eq finished CLEAN (T=298.5, density 0.994, P~0) -- prod died step-1 with
#            "Periodic box dimensions changed too much" on all 5 restarts. Root cause:
#            (A) boresch.RST reference r0=5.135A but actual L1-P1 ~65A in start coords
#            (RESTRAINT=50,022 kcal/mol at lambda0) -> giant force -> GPU step-1 abort.
#            New _verify_boresch_reference() re-measures L1-P1 vs the SAME start frame
#            (system.inpcrd/min.rst) with minimum image and FAILS LOUDLY at write time
#            if |actual-r0|>2A or > half-box (periodic-image anchor). (B) prod.in had
#            vlimit=-1 (OFF) -> added vlimit via _vlimit_block(). (C) box-drift retry
#            re-ran prod from the SAME eq.rst 5x (futile); now inserts ONE short GPU NPT
#            regrid step (dt=0.001, vlimit, 50k steps) that writes a fresh restart and
#            resumes prod -- stays on GPU. (D) boresch-gate now HARD-STOPS (exit 3) when
#            amber_md is not importable on the node, instead of running UNVALIDATED.
#   2.5.31c- ROOT-CAUSE FIX: heat &wt keyword 'TEMP_0' must be 'TEMP0'. pmemd aborted
#            EVERY leg with "Invalid TYPE flag" at heat (proven from real heat.out:
#            restraint/decharge/vdw all identical). Prior 2.5.29/2.5.31 chased the
#            &wt FORMAT; the defect was the spelling (manual 18.6.1: nmropt=1 TEMP0
#            IS supported under TI). Fixed generator fep._heat_in, FEP-SPell template,
#            sample + smoke-test. Validator: new check_wt_type_keywords() makes any
#            non-canonical &wt type= a FATAL gate error (TEMP_0->TEMP0 hint). Full
#            tree audited: no other non-canonical &wt cards.
#   2.5.31b- VALIDATOR TIMING/NOISE FIX + VERSION-TRACK DISCIPLINE.
#            (a) Two-phase driver gate: content checks pre-script, then the -ref
#                run-script check AFTER run_<leg>.lsf is written -- both before any
#                bsub. Fixes the 2.5.31a gap where the critical missing-ref check
#                was inert (script not yet written) and warnings spammed 16-44x.
#            (b) check_run_script runs ONCE per leg, targets run_<leg>.lsf, ignores
#                analyze_*/cycle_close_* scripts.
#            (c) VERSION TRACK now updated EVERYWHERE on each bump: __init__.py,
#                VERSION, run_amber.py, README.md banner, and the GUI Home.py
#                "What's new" page (which had been stale since 2.5.26). Added
#                tools/check_version_sync.py (CI-style guard) and tools/bump_version.py.
#                NOTE re the lig_12944901 crash: the .lsf was stamped 2.5.31a but the
#                heat.in carried the pre-2.5.29 inline &wt -> the COMPUTE ENV imported
#                a stale amber_md. Reinstall this package into that env so the
#                IMPORTED code matches the stamp (see README "Install / version sync").
#   2.5.31a- VALIDATOR FALSE-POSITIVE FIX (found by the gate on the first real ABFE
#            run): (1) the inline "&wt type='END' /" terminator written by
#            _restraint_block in min/dens/eq/prod is TOLERATED by pmemd (years of
#            working decharge/vdw legs) -- only DATA-bearing inline &wt cards (e.g.
#            TEMP_0) are fatal; stop flagging the bare END terminator. (2) timask2
#            is legitimately EMPTY for ABFE (single-region decoupling: the driver
#            emits timask2=''); require timask2 only when scmask2 is non-empty
#            (two-region/RBFE). The gate now PASSES a valid ABFE leg while still
#            blocking the real inline-TEMP_0 and missing-timask2(RBFE) cases.
#   2.5.31 - WIRED-IN VALIDATION: mdin checks are now an importable module
#            amber_md/mdin_validator.py AND a hard GATE in fep_driver.run_fep --
#            after each leg's .in files + run script are written and BEFORE any
#            bsub. A fatal issue aborts the run (rc=2) and logs the cause; nothing
#            malformed reaches the GPU. tools/validate_mdin.py imports the shared
#            module. New --skip-mdin-check overrides (not recommended). GUI shows a
#            note for Amber ABFE/RBFE.
#   2.5.30 - TOOLING: added tools/validate_mdin.py, a GPU-free pre-submit validator
#            that parses every generated heat/dens/eq/prod .in against the pmemd
#            rules that have repeatedly aborted runs at SETUP (one defect at a time):
#              * inline "&wt type=..." weight card  -> "Invalid TYPE flag"  (2.5.29)
#              * "*" residue wildcard in restraintmask -> "unknown symbol:*" (2.5.28)
#              * ntr=1 stage whose run command lacks -ref -> "refc" open err (2.5.27)
#            plus heat-is-NVT and DISANG-file-exists checks. Verified to FAIL a
#            known-bad window (all 4 issues) and PASS a clean v2.5.29 window.
#            Run: python tools/validate_mdin.py <leg_or_run_dir>  (exit!=0 blocks).
#   2.5.29 - HOTFIX (regression, ALL legs): the v2.5.25 heat stage emitted the &wt
#            temperature-ramp card inline ("&wt type='TEMP_0', ...") which pmemd's
#            weight-change reader rejects -> "Error: Invalid TYPE flag in line:
#            TEMP_0 ..." at the heat stage of EVERY leg (decharge/vdw legs that
#            used to run failed too). Reformatted the heat &wt block to the proven
#            FEP-SPell-ABFE layout (equilibration_fepspell.py lines 86-93): "&wt"
#            alone on its line, then "type = 'TEMP_0',". The dens/eq/prod END-only
#            &wt tail (_restraint_block) is unchanged -- it has run for years.
#   2.5.28 - HOTFIX (restraint mask + heat template): pmemd group-input parser
#            (ntr=1 restraintmask) rejects '*' wildcards -> heat aborted with
#            "unknown symbol:*". Replaced Cl*/Na* with explicit ion names; verified
#            via parmed/AmberMask vs the real 272,652-atom complex_restraint.parm7
#            (same 6,850 atoms; :1/:2=LIG/LIG 37 atoms matched; Boresch on LIG+prot).
#            Also rewrote _heat_in to emit a single nmropt=1 and single &wt END.
#   2.5.27 - HOTFIX (regression from 2.5.25): restrained stages (heat/dens/eq)
#            emit ntr=1 but the run_stage commands never passed -ref, so pmemd
#            aborted at "5. REFERENCE ATOM COORDINATES" with
#            "Unit 10 Error on OPEN: refc" (seen first on complex_restraint heat,
#            lambda=0.000). Added -ref <stage-start-coords> to all six restrained
#            run_stage lines (heat/dens/eq x array+HREMD). min stays ntr=0/no-ref.
#   2.5.26 - CLEAN-RUN READINESS: heat-stage stability gate (peak T >
#            heat_temp_max_K -> exit 70 BEFORE dens) in both LSF-array and
#            HREMD run paths, mirroring the eq gate; tools/verify_stage_order.sh
#            GPU-free pre-submit dry-checker; README_CLEAN_RUN.md.
#   2.5.25 - RELIABILITY (decouple/vdw legs): added explicit NVT HEAT stage
#            (ntb=1,ntp=0,dt=0.001,TEMP_0 ramp,restrained) between min and dens
#            so the barostat never acts while heating a half-decoupled soft-core
#            ligand (root cause of mid-lambda T->15,000 K blow-ups); default
#            positional restraint on dens+eq (posres_mask_default); vlimit valve
#            on heat/dens/eq. New stage order: min->heat->dens->eq->prod.
#            Cross-checked vs FEP-SPell-ABFE/BAT.py/GHOAT.py; soft-core block
#            already matched FEP-SPell-ABFE exactly. See
#            CHANGES_v2.5.25_decouple_vdw_reliability.md.
#   2.5.24 - FIX (restraint leg): octahedral ChBox restore (solvateOct+'set box'
#            reset angles to 90 deg -> pmemd NPT crash); build dual-copy topology
#            from EQUILIBRATED coords; inherit neutral system (no re-ionization);
#            drop crgmask on the restraint leg (pure :1/:2 TI, ifsc=0) to match
#            the validated FEP-SPell-ABFE leg. Independent rk->2*rk + no-crgmask
#            confirmation from BAT.py. See
#            CHANGES_v2.5.24_restraint_chbox_and_reference_alignment.md.
#   2.5.23 - FIX (restraint leg, Option B, ARCHITECTURE): the four prior
#            single-prmtop attempts (2.5.20-2.5.22) all failed because the
#            complex_restraint leg reused complex.prmtop (ligand=:872, last
#            residue) and tried to ramp ONLY a Boresch potential with both TI
#            ends sharing the same fully-interacting ligand. Now a DEDICATED
#            dual-copy topology (combine { LIG LIG protein }: :1 real, :2
#            dummy) runs ORDINARY TI -- timask1=:1/timask2=:2, ifsc=0,
#            crgmask=:2, fixed-k Boresch. New module abfe_restraint_topology
#            builds it + remaps Boresch atom indices. See
#            CHANGES_v2.5.23_restraint_dualcopy.md.
#   2.5.22 - BUGFIX (restraint leg, take 3, physics): ifsc=0 made the
#            fully-interacting ligand a HARD-core TI region -> VDWAALS/EEL
#            overflow at step 1 (same structure decharged fine). Restraint
#            leg now runs ifsc=1 + GTI soft-core with identical :LIG TI/sc
#            end states (zero molecular dV/dl); restraint lambda-scaling is
#            unchanged so the leg still integrates the restraint work.
#   2.5.21 - BUGFIX (restraint leg, take 2): empty TI masks tripped pmemd
#            "must match at least one atom". Restraint leg now uses identical
#            non-empty TI masks (ligand on both ends, zero molecular dV/dl)
#            AND lambda-scales the Boresch force constants per window so the
#            leg actually integrates the restraint work.
#   2.5.20 - BUGFIX: complex_restraint leg crashed at min (timask1/2 atom
#            count mismatch with ifsc=0). Restraint leg now uses empty,
#            matching TI masks; the Boresch potential is applied via
#            nmropt/DISANG, not a topology perturbation.
#   2.5.19 - Cleanup + equilibration hardening: removed dead fep_mdin_v2516;
#            positional restraints carried through dens+eq (opt-in posres_mask).
#   2.5.18 - Wired the FEP-SPell ABFE restraint/correction/charge-correction
#            into the PRODUCTION path (fep.py/fep_driver). Boresch restraint now
#            self-consistent (6/6 DOF); Deng-Roux correction; auto charge corr.
#   2.5.17 - ABFE restraint self-consistency + automation (FEP-SPell port).
#            Corrected Boresch atom ordering; NEW measure/auto-select/staged-eq;
#            post-equilibration fail-fast + self-repair gate
#            (abfe_integration.verify_or_reselect_boresch); manual masks
#            override auto-selection. See CHANGES_abfe_v2.5.17.md.
#   2.5.16 - ABFE protocol + free-energy bias fixes.
#   2.5.15 - FIX (final71): WINDOW RECOVERY + MBAR rank guard + dense vdw.
#            (A) RECOVERY (self-heal): detonated windows could not recover --
#                the rerun resumed PRODUCTION from the blown-up eq.rst, so it
#                re-detonated every attempt. Fixes:
#                  * recovery_stages(): for a PHYSICS blow-up (blowup_temperature,
#                    shake_failure, box_drift, softcore_com, unknown,
#                    missing_no_error) rebuild dens->eq->prod from a CLEAN coord
#                    source (min.rst / system.inpcrd) so the remediation edits to
#                    dens.in/eq.in/prod.in are actually exercised.
#                  * purge_poisoned_restarts(): delete eq.rst/dens.rst/prod.* (or
#                    just prod.* for transient/external kills) before rerun so
#                    pmemd cannot resume (irest=1) from detonated coordinates.
#                  * GUARANTEED last-resort: on the final attempt apply a very
#                    gentle fixed protocol (dt=0.0005, gamma_ln=10, Berendsen
#                    barostat, softened soft-core scalpha/scbeta) from clean
#                    coords so a window completes rather than looping.
#                  * external_kill keeps the cheap resume-from-good-eq.rst path.
#            (B) RANK GUARD (analysis): SVD rank/cond of u_nk; MBAR only fit when
#                rank==n_states and cond<=1e8, else flagged ill-posed (the -1477
#                artifact). BAR/TI unaffected.
#            (C) vdw_lambdas_dense: opt-in 40-window schedule (use_dense_vdw) for
#                better overlap/conditioning. Legacy 28-window kept as default.
#   2.5.14 - FIX (final70): root-caused MBAR failure on REAL full-leg data.
#            (1) RAGGED MBAR GRID killer: a vdw leg built from TWO lambda
#                schedules (23 old + 5 refined = 28-grid) -> non-square u_nk
#                that MBAR/BAR cannot solve; unfixable post-hoc. Analyzer now
#                DETECTS inconsistent per-window MBAR grids and ABORTS the leg
#                with a clear message instead of emitting silent garbage.
#            (2) MIScalibrated sanitizer: _UNK_SANE_MAX_KT=1e4 + column/row
#                DROP decimated CLEAN matrices (~95% frames lost, 5x worse dG).
#                Real corruption is only the decoupled end-state singularity
#                (non-finite + ~1e9); Boltzmann-correct fix is CLIP-to-ceiling
#                (zero weight), never drop. Recovers truth to 0.03 kcal/mol
#                under 22% corruption. _UNK_SANE_MAX_KT=1e3 (reduced units).
#   2.5.13 - FIX: MBAR/BAR still diverged (dG=-710 vs TI=+8.9) on COMPLETE legs
#            after 2.5.12. Real-data root cause (lig solvent_vdw/l=0.900
#            prod.out): the soft-core decoupled END STATE (e.g. the lambda=0.0
#            column) overflows to a Fortran "****" field in MOST frames, while
#            every OTHER state column is clean. The 2.5.12 "drop the whole
#            frame if ANY cell is bad" guard therefore dropped ~all frames
#            (observed 32340), emptied half the lambda columns, and pymbar hit
#            DLASCL/SVD on the degenerate matrix. Fix: COLUMN-AWARE sanitation
#            -- drop only the unrecoverable state column(s) (bad in >50% of
#            frames), keep all good frames/states; never median-fill (that
#            fabricates samples and biases dG by ~3.7 kcal/mol). Validated to
#            0.0000 kcal/mol vs ground truth on a synthetic harmonic ladder
#            with a 90%-corrupted end column. TI unaffected. Analysis-only
#            change: re-run with --resume --analyze; no MD re-run required for
#            legs whose windows already completed.
#   2.5.12 - FIX: MBAR/BAR returned garbage (dG = -100..-1150 kcal/mol) on
#            EVERY leg, including fully-complete ones, while TI was correct.
#            Root cause: soft-core end-state overflow blows up individual u_nk
#            cells to ~1e6 kT; the final60 clamp-to-row-median guard failed to
#            bound them (an all-overflow row has a huge median), so pymbar got
#            a corrupt reduced-potential matrix -> DLASCL/SVD failure. Fix:
#            drop overflow frames on an absolute kT bound, per-lambda, keeping
#            the matrix square. TI unaffected (uses dHdl).
#   2.5.11 - Refined vdW lambda schedule (config default). Under uniform 0.05
#            spacing, complex_vdw windows at lambda 0.70/0.75 hit eq
#            instability (exit 71) and several mid-lambda windows hit prod box
#            drift (exit 255) -- the soft-core danger zone. Halved the spacing
#            across 0.6-0.85 (added 0.625/0.675/0.725/0.775/0.825; 23->28
#            windows). NOT a code bug -- the eq gate / box-drift handling did
#            their job; this reduces the per-window perturbation so they don't
#            trigger. Changing the schedule requires a full vdw re-run.
#   2.5.10 - FIX: prod box-drift auto-restart never triggered. pmemd.cuda writes
#            'Periodic box dimensions have changed too much...' to STDOUT, but
#            run_prod_with_restart grepped prod.out (the -o mdout), which never
#            contains it -> a box-drift rc=255 was misclassified as a non-box
#            error and the window failed on attempt 1 with no retry. Now pmemd
#            stdout+stderr is tee'd to prod.console.<attempt> and the drift
#            detection greps that AND prod.out. Verified the retry loop engages
#            and recovers. (array builder only -- HREMD has no per-window prod
#            restart by design; bash -n release gate still passes.)
#   2.5.9  - CRITICAL FIX: the v2.5.7 eq-gate patch left a DUPLICATED
#            'if [ -n "$maxT" ]; then' line in build_lsf_array, so the
#            generated .lsf had an unbalanced if/fi -> bash 'syntax error:
#            unexpected end of file' -> every window exited 2 after eq, before
#            prod. Removed the duplicate. The generator output is now validated
#            with 'bash -n' (syntax) for both array + HREMD builders so a
#            malformed script can never ship again.
#   2.5.8  - True host-avoidance for genuine node failures. HPCConfig.avoid_hosts
#            is a persistent GPU blocklist emitted as -R "select[hname!=...]"
#            into every generated GPU job (initial array + HREMD launch). The
#            self-heal bsub rerun additionally excludes any host that already
#            ran/killed the window (parsed from the LSF .out), plus a
#            --avoid-hosts CLI blocklist; --no-host-avoidance opts out. Also
#            fixed the self-heal 'bsub -K < script' stdin redirect.
#   2.5.7  - CRITICAL FIX: the final60 eq temperature-stability gate computed
#            peak T via 'grep ... | sort -gr | head -n1'. Under 'set -o
#            pipefail' + 'set -e', head's early pipe-close made upstream grep/
#            sort die with SIGPIPE (rc=141), aborting EVERY window on EVERY
#            node right after eq, before prod -> no prod.out. Now computed in a
#            pipefail/errexit-free subshell with a single-pass awk max. Fixes
#            the all-windows exit-141 failure. (array + HREMD builders.)
#   2.5.6  - Version stamped into every generated .lsf header and shown on the
#            GUI homepage; amber_md.version helpers. self-heal classifies
#            external scheduler/node kills as `external_kill` -> clean resubmit.
#   2.5.5  - final61: per-window self-heal (abfe_self_heal[/ _cli]) wired into
#            fep_driver via --self-heal; final60 eq-stability gate + prod
#            box-drift auto-restart; FEP estimator sanitizer/reconcile.
#   2.5.0  - Unified Setup & Launch wizard; OpenMM MM-GBSA; atom-mapping
#            inspection; promote MM-GBSA -> FEP; config round-trip fix.
# ---------------------------------------------------------------------------
__version__ = "2.6.0"

# v2.5.31g: PEP 440-legal projection of __version__ for packaging tools.
def _pep440(v):
    import re as _re
    m = _re.match(r"^(\d+(?:\.\d+)*)([a-z])?$", v.strip())
    if not m:
        return v
    base, letter = m.group(1), m.group(2)
    if not letter:
        return base
    return f"{base}.post{ord(letter) - ord('a') + 1}"

__pep440_version__ = _pep440(__version__)

# Short build/codename shown alongside the semver (helps disambiguate the
# directory name from the code version when reading a generated .lsf).
__build__ = "gui-streamlined"
