"""Alchemical FEP / TI (dual-topology softcore) for Amber.

v2.4.23 (patch set):
  * FIX: decharge leg now emits scalpha/scbeta whenever ifsc>0.
    ABFE single-topology forces ifsc=1 on the decharge stage, but
    _stage_softcore_block() previously returned "" for decharge, so
    pmemd used scalpha=scbeta=0 -> a HARD soft-core core. Benign in
    solvent (free ligand), but in the crowded complex it produced a
    dV/dlambda singularity (dV/dl -> -113 at high lambda, MBAR spread
    ~2700, dG_bind unphysical +63.5). Softening is now applied to both
    decharge and vdw stages whenever softcore is active.
  * FEATURE: analyzer CSVs are now self-identifying. dvdl_summary.csv
    and dG_estimators.csv gain a leading "# leg: <name>" comment AND a
    `leg` column (e.g. complex_decharge), so solvent vs complex files
    can no longer be confused.

v2.4.22 (patch set):
  * FIX: dedicated density pre-equilibration stage (dens) inserted
    between min and eq (NPT, taup=1.0); eq now continues from dens.rst.
    Converges box volume before production so pmemd.cuda no longer
    halts with "Periodic box dimensions have changed too much" (rc=255).
  * FIX: prod runs via run_prod_with_restart() -- on the GPU box-drift
    error it auto-restarts from the latest restart file (up to 5x),
    Amber's prescribed recovery; any other failure surfaces normally.

v2.4.21 (patch set):
  * FIX: TI-anchored headline selection in FEPAnalyzer.run() --
    BAR/MBAR are used as headline ONLY if they agree with TI within
    _ESTIMATOR_CONSISTENCY_KCAL (2.0). Stops broken MBAR (e.g. -130,
    -150 on the solvent legs) from winning the headline just because
    |dG| < 500. Records estimator_spread_kcal / estimator_converged.
  * FIX: cycle-closer trust gate now requires (a) all legs complete,
    (b) every leg estimator_converged, and (c) |dG_bind| <=
    _DG_BIND_SANITY_KCAL (25). A complete-but-unphysical result is no
    longer stamped trusted=true. Reports the precise reason.

v2.4.20 (patch set):
  * FIX: add top-level `import re` (build_lsf_array_resume used
    re.sub but re was never imported -> NameError on --resume).
  * FIX (solvent vdw): _prod_in / _eq_in emit nscm=0. pmemd.cuda
    aborted at prod init (rc=255, no NaN) for FULLY soft-core windows
    ("Molecule N is completely softcore and skipped for C.O.M.")
    because COM-motion removal trips on a fully decoupled molecule
    restarting with velocities (irest=1,ntx=5). nscm=0 disables that
    step so high-lambda decoupled vdw windows restart cleanly.
  * FIX (complex vdw / the cluster): build_lsf_array + build_lsf_hremd emit
    span[hosts=1] and rusage[mem=fep_mem_mb] (default 8192 MB),
    matching LSFSubmitter._header. Prevents host-OOM that killed eq
    mid-run (truncated, no eq.rst) and complies with mandatory Indy
    the cluster memory-request enforcement.
  * FIX: _eq_in writes ntwr (intermediate restart) so a killed eq
    leaves an eq.rst to resume from instead of restarting from scratch.

v2.4.19 (patch set):
  * FIX Bug 1: submit_cycle_closer() reads two-stage analyzer JID keys
    -> cycle-closer is actually submitted and ABFE_RESULT.txt written.
  * FIX Bug 3: build_cycle_closer_lsf() reads decharge/vdw schedules via
    getattr(..., cfg.lambdas) so a sparse config cannot crash the build.
  * FIX Bug 4: cycle-closer _run() only trusts a COMPLETE cached
    summary.json; otherwise re-analyzes (resume-safe).

v2.4.15:
  * FIX: decharge stage (ifsc=0) no longer emits scmask1/scmask2.
    pmemd aborts with "scmask1 is only used with ifsc > 0!" when a
    softcore mask is present while ifsc=0, which killed every window
    of both decharge legs. scmask is now emitted only when ifsc>0,
    via the new FEPSetup._scmask_block() helper, mirroring how
    _crgmask_block() / _stage_softcore_block() are already gated.
  * FIX: FEPAnalyzer no longer reports a clean dG from an INCOMPLETE
    leg. It tracks n_requested vs n_windows, records
    missing_windows/complete, refuses to mark an incomplete leg's
    estimator as the trusted "*" headline, and writes a `complete`
    column + a leading "# INCOMPLETE ..." comment into the CSVs.
  * FIX: _parse_dvdl hardened - scans ALL DV/DL lines in the averages
    block, tolerates "="-less formats, and rejects non-finite values.
  * FIX: build_lsf_array() now runs each pmemd stage through a
    run_stage() shell wrapper that, on failure, greps the stage mdout
    for the real cause (STOP PMEMD / NaN / vlimit exceeded / input
    errors) and prints a "STAGE FAILED" diagnostic before exiting with
    that stage's code -- instead of set -e dying silently as exit 1.
  * FIX: build_lsf_hremd() gained the same run_stage() diagnostics for
    the per-window min/eq prep and a guarded mpirun replica-exchange
    call that, on failure, dumps rem.log + each window prod.out tail.
  * FIX: analyzer LSF (set -uo pipefail, no -e) now propagates the
    Python exit status and prints an explicit OK/INCOMPLETE/FAILED
    banner, so a crashed analyzer no longer reports success to the
    dependent cycle-closer.
  * FIX: cycle-closer prints per-leg completeness, refuses to emit a
    trusted dG_bind when any leg is missing/incomplete, and exits
    nonzero so its LSF record reflects the real outcome.
  * FIX: submit_cycle_closer() now depends on ended(analyzer) instead
    of done(analyzer). Since analyzers exit nonzero on INCOMPLETE/
    FAILED legs, done() would have blocked the cycle-closer from ever
    launching; ended() lets it always run and write a definitive
    (possibly UNTRUSTED) ABFE_RESULT instead of silently never firing.

v2.4.16 (alchemical-scheme correctness):
  * FIX: ABFE decharge stage now runs ifsc=1 (was ifsc=0). With the
    single-topology mask (timask1=:LIG, timask2='') the atom counts
    are intentionally mismatched, which is ONLY legal under softcore;
    ifsc=0 made pmemd abort "timask1/2 must match the same number of
    atoms for non-softcore run". Decharging is driven by crgmask, not
    by turning softcore off.
  * FIX: crgmask is now emitted on the DECHARGE stage (it zeroes the
    ligand charges) instead of the vdw stage, where it was both wrong
    and redundant.
  * FIX: _mask_block() now produces the ABFE single-topology masks
    (timask2='', scmask stage-aware) for EVERY leg when an ABFE/
    two-stage decoupling is in effect -- previously it keyed off
    self.boresch, so the solvent leg (no Boresch restraint) fell back
    to RBFE defaults and emitted timask2=':MOD', a residue that does
    not exist in the solvent topology (matched 0 atoms).

v2.4.18 (resume support):
  * NEW: incomplete_indices() / submit_leg_resume() -- detect windows
    whose prod.out lacks pmemd's completion marker (killed by HPC
    maintenance/preemption) and resubmit ONLY those as a sparse LSF
    array, then re-wire the analyzer so the cycle-closer fires and
    ABFE_RESULT.txt is written. Finished windows are never rerun.
    submit_leg_resume() mirrors submit_leg()'s exact _bsub_submit
    usage (walltime required; dependency via extra_args=[-w ended()]).
    Drives the GUI "Resume" button.

v2.4.17 (decharge/vdw scheme corrected to single-prmtop softcore):
  * CONTEXT: abfe_topology.py builds ONE ordinary charged prmtop and
    fep_driver.py feeds that SAME prmtop to both the decharge and vdw
    legs (no dual-topology / no separately-built decharged prmtop).
    That is exactly Amber manual sec.27.1.8.2 "Absolute free energy
    using soft core": icfe=1, ifsc=1, timask1=:LIG, scmask1=:LIG,
    timask2='', scmask2='' -- clambda 0->1 decouples the ligand.
  * FIX: v2.4.16 moved crgmask to the DECHARGE stage and emptied its
    scmask. That is the DUAL-topology idiom and is WRONG here: with a
    single charged prmtop, a decharge-stage crgmask zeroes :LIG in
    BOTH end states, making V0==V1 so the stage integrates to ~0 (a
    silent, plausible-but-wrong result). Reverted to the correct
    single-prmtop stepwise scheme:
      - decharge: ifsc=1, scmask1=:LIG, NO crgmask  -> clambda removes
        the ligand CHARGES via softcore electrostatics (manual Eq.27.7)
      - vdw:      ifsc=1, scmask1=:LIG, crgmask=:LIG -> charges held at
        zero (manual sec.27.1.8.3 note: "to set up a soft core vdW
        transformation, the flag crgmask can be added") while clambda
        removes the ligand vdW.
  * KEPT from v2.4.16: ifsc=1 on both stages (the ifsc=0 crash fix);
    single-topology timask2='' on every leg (the :MOD fix); solvent
    leg decouples the same :LIG as the complex leg.

v2.4.14:
  * NEW: FEPSetup.build_cycle_closer_lsf() — single-CPU LSF job that
    waits on BOTH per-leg analyzers, then writes fep/ABFE_RESULT.txt +
    fep/ABFE_RESULT.json via ΔG_bind = -(ΔG_complex+restr - ΔG_solvent).
    The full ABFE pipeline now completes with no Python or human alive
    after submit.
  * NEW: FEPSetup.submit_cycle_closer() — wraps the bsub call.

v2.4.12 PATCHES (kept):
  A. submit_leg() also submits a dependent analyzer LSF job.
  B. FEPAnalyzer.run() rejects insane MBAR (|dG|>500 or NaN err),
     records headline_estimator in the result.
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
import math
import re
from .submit import _bsub_submit
from .utils import run, ensure_dir
from .logger import get_logger
from .config import HPCConfig, MDConfig, FEPConfig
from .version import lsf_banner, runtime_echo  # v2.5.6: stamp version into .lsf
log = get_logger()

_RESTRAINT_HEADER = (
    "# Boresch-style 6 restraints (1 distance, 2 angles, 3 dihedrals)\n"
    "# Generated by amber_md.fep.FEPSetup\n"
)
_HEADLINE_SANITY_MAX_KCAL = 500.0
_ESTIMATOR_CONSISTENCY_KCAL = 2.0   # v2.4.21: BAR/MBAR must agree with TI within this
_UNK_SANE_MAX_KT = 1.0e3        # v2.5.14: u_nk REDUCED-energy ceiling (kT). Values reaching
                               # the sanitizer are alchemlyb RELATIVE reduced potentials
                               # (physical median ~O(1-100) kT). End-state singularities
                               # print as **** (non-finite) or finite ~1e9 garbage. 1e3
                               # sits above any physical |u_nk| yet below the overflow
                               # tail, so genuine cells are untouched and ONLY true
                               # singularities are clipped. (Old 1e4 ABSOLUTE cap +
                               # col/row DROP decimated clean matrices -- removed.)
_UNK_COND_MAX = 1.0e8           # v2.5.15: max u_nk cond before MBAR ill-posed.
_UNK_RANK_RTOL = 1.0e-6        # v2.5.15: rel SVD tol for u_nk numerical rank.
_DVDL_SANE_MAX_KCAL = 5.0e3     # final39: TI dV/dl endpoint-singularity cap (kcal/mol).
_DVDL_OUTLIER_MAD = 25.0        # final39: per-lambda robust outlier cut (x MAD) on dV/dl frames
_ERROR_BLOCKS = 5               # v2.4.25 block-average error (BAT.py-style)



class FEPSetup:
    def __init__(self, work_dir, cfg, md_cfg, hpc_cfg,
                 hremd=False, exchange_freq=1000, boresch=None):
        self.root          = ensure_dir(Path(work_dir) / "fep")
        self.cfg           = cfg
        self.md_cfg        = md_cfg
        self.hpc_cfg       = hpc_cfg
        self.hremd         = bool(hremd)
        self.exchange_freq = int(exchange_freq)
        self.boresch       = boresch
        self._stage        = None
        self._suppress_correction = False

    def _restraint_block(self):
        if not self.boresch: return "", ""
        return "  nmropt=1,\n", "&wt type='END' /\nDISANG=boresch.RST\n"

    def _is_abfe(self):
        """True for an absolute (single-ligand decoupling) calculation.

        ABFE is signalled by the two-stage decharge/vdw protocol (the
        presence of crgmask / two_stage), NOT by Boresch restraints:
        the solvent leg of an ABFE has no Boresch restraint but is still
        absolute. Keying off boresch (the old behaviour) made the
        solvent leg fall through to RBFE defaults (timask2=':MOD')."""
        if getattr(self.cfg, "two_stage", False):
            return True
        if self._stage in ("decharge", "vdw"):
            return True
        return bool(self.boresch)

    def _ligand_mask(self):
        """The single perturbed-ligand mask for ABFE (e.g. ':LIG')."""
        return getattr(self.cfg, "crgmask", None) or self.cfg.timask1

    def _mask_block(self):
        """Return (timask1, timask2, scmask1, scmask2).

        ABFE  -> single-topology decoupling of one ligand (manual
                 sec.27.1.8.2). The SAME charged prmtop is used for both
                 stages, so the softcore region is :LIG in BOTH:
                   timask1=:LIG, timask2='', scmask1=:LIG, scmask2=''
                 The stages differ only by crgmask (see _crgmask_block):
                   decharge -> no crgmask (clambda removes charges)
                   vdw      -> crgmask=:LIG (charges held off, clambda
                               removes vdW)
        RBFE  -> dual-topology A->B from config masks (unchanged).
        """
        # v2.5.23 FIX (Option B dual-copy): the restraint leg now runs on a
        # DEDICATED topology built by abfe_restraint_topology.build_restraint_
        # topology(), laid out as  combine { LIG LIG protein }:
        #     :1 = real ligand   (TI end state V0)
        #     :2 = dummy ligand  (TI end state V1, charges zeroed via crgmask)
        # This is ordinary TI with MATCHING 37-atom TI regions, ifsc=0, and a
        # FIXED-k Boresch restraint (no lambda scaling). It replaces the four
        # failed single-prmtop attempts (v2.5.20 empty-mask atom-count abort,
        # v2.5.21 ifsc=0 hard-core VDWAALS overflow, v2.5.22 ifsc=1 soft-core
        # SC_VDW overflow), all of which tried to ramp only a restraint while
        # both TI ends shared the same fully-interacting ligand. See
        # abfe_restraint_topology.py for the topology + Boresch index remap.
        if self._stage == "restraint":
            # v2.5.33: dual-copy TI (FEP-SPell-faithful). :1 = real ligand,
            # :2 = identical second copy; clambda transforms :1<->:2 with FIXED-k
            # Boresch restraining :1. ifsc=0 (parameter-identical copies).
            t1 = getattr(self.cfg, "restraint_timask1", ":1")
            t2 = getattr(self.cfg, "restraint_timask2", ":2")
            return (t1, t2, "", "")
        if self._is_abfe():
            lig = self._ligand_mask()
            return (lig, "", lig, "")
        return (self.cfg.timask1, self.cfg.timask2,
                self.cfg.scmask1, self.cfg.scmask2)

    def _softcore_block(self):
        """Complete GTI soft-core + smoothing block (v2.5.16).

        Matches the published FEP-SPell-ABFE protocol. The previous block
        emitted only scalpha/scbeta + gti_add_sc/gti_chg_keep, leaving
        gti_ele_sc / gti_vdw_sc / gti_lam_sch UNSET -> the half-configured
        TI soft-core that diverged at mid-lambda (T~15000 K, box drift).
        """
        c = self.cfg
        g = lambda n, d: getattr(c, n, d)
        # v2.5.70: STAGE-AWARE vdW soft-core + charge handling. On the DECHARGE
        # stage the ligand must decharge with HARD vdW (gti_vdw_sc=0) and charges
        # genuinely removed (gti_chg_keep=0); otherwise opposite partial charges
        # approach through a softened core -> charging singularity (dV/dl ~ -507
        # kT, **** in u_nk, MBAR SVD failure). The VDW stage keeps full soft-core
        # (charges already held off by crgmask=:LIG) -- so it is unchanged.
        if self._stage == "decharge":
            _vdw_sc = g('decharge_gti_vdw_sc', 0)
            _chg_keep = g('decharge_gti_chg_keep', 0)
        else:
            _vdw_sc = g('gti_vdw_sc', 1)
            _chg_keep = g('gti_chg_keep', 1)
        return (
            f"  scalpha={g('scalpha',0.2)}, scbeta={g('scbeta',50.0)},\n"
            f"  gti_lam_sch={g('gti_lam_sch',1)}, "
            f"gti_ele_sc={g('gti_ele_sc',1)}, gti_vdw_sc={_vdw_sc},\n"
            f"  gti_scale_beta={g('gti_scale_beta',0)}, "
            f"gti_cut_sc={g('gti_cut_sc',0)}, gti_cut={g('gti_cut',1)},\n"
            f"  gti_add_sc={g('gti_add_sc',1)}, "
            f"gti_chg_keep={_chg_keep}, "
            f"tishake={g('tishake',1)}, logdvdl={g('logdvdl',0)},\n"
        )

    def _icfe_value(self):
        """v2.5.36 (Option A): restraint leg ramps ONLY the Boresch potential
        via lambda-scaled &rst force constants (no alchemical TI) -> icfe=0."""
        return 0 if self._stage == "restraint" else 1

    def _ti_kw_line(self, cl):
        if self._stage == "restraint":
            return ""
        t1, t2, s1, s2 = self._mask_block()
        return (f"  icfe=1, ifsc={self._ifsc_value()}, clambda={cl},\n"
                f"  timask1='{t1}', timask2='{t2}',\n"
                f"{self._scmask_block()}{self._crgmask_block()}{self._stage_softcore_block()}")

    def _ti_kw_line_mbar(self, cl):
        if self._stage == "restraint":
            return ""
        t1, t2, s1, s2 = self._mask_block()
        return (f"  icfe=1, ifsc={self._ifsc_value()}, clambda={cl}, ifmbar=1,\n"
                f"  mbar_states={len(self._active_lambdas)},\n"
                f"  mbar_lambda={','.join(f'{l:.4f}' for l in self._active_lambdas)},\n"
                f"  timask1='{t1}', timask2='{t2}',\n"
                f"{self._scmask_block()}{self._crgmask_block()}{self._stage_softcore_block()}")

    def _ifsc_value(self):
        # ABFE single-topology masks have mismatched atom counts
        # (timask1=:LIG vs timask2=''), which is legal ONLY with
        # softcore enabled -- so BOTH decharge and vdw use ifsc=1.
        # (RBFE matched-atom dual topology may use ifsc=0.)
        # v2.5.16: the restraint-removal leg ramps ONLY the Boresch potential
        # on the fully-interacting complex; the ligand keeps full vdW+charge,
        # so there is NO soft-core region -> ifsc=0.
        if self._stage == "restraint":
            # v2.5.33: dual-copy :1/:2 TI -- parameter-identical copies, no soft
            # core needed -> ifsc=0 (FEP-SPell setting).
            return 0
        if self._is_abfe():
            return 1
        return 0 if self._stage == "decharge" else 1

    def _stage_softcore_block(self):
        # v2.4.23: scalpha/scbeta must be present for BOTH stages whenever
        # ifsc=1. ABFE single-topology forces ifsc=1 on the decharge leg
        # (mismatched timask atom counts), so the decharge leg ALSO needs
        # softening -- otherwise scalpha=scbeta=0 (pmemd default) gives a
        # hard core that blows up dV/dlambda in the crowded complex
        # (dV/dl -> -113, MBAR spread 2700). Benign in solvent (free
        # ligand rarely overlaps), so this is safe for all legs.
        if self._ifsc_value() > 0:
            return self._softcore_block()
        return ""

    def _scmask_block(self):
        """Emit scmask only when (a) softcore is on AND (b) a non-empty
        softcore region is defined. pmemd aborts with
        'scmask1 is only used with ifsc > 0!' if scmask is present while
        ifsc=0, and an ABFE decharge stage deliberately has an EMPTY
        scmask (charges-only), so the line must be omitted there too."""
        if self._ifsc_value() <= 0:
            return ""
        _, _, s1, s2 = self._mask_block()
        if not s1 and not s2:
            return ""
        return f"  scmask1='{s1}', scmask2='{s2}',\n"

    def _crgmask_block(self):
        # Single charged prmtop, stepwise protocol (manual sec.27.1.8.3
        # note): the DECHARGE stage removes :LIG charges via clambda and
        # must NOT pin them with crgmask. The VDW stage keeps the charges
        # held at zero (crgmask=:LIG) while clambda removes vdW.
        # v2.5.16: the restraint leg perturbs neither charge nor vdW (only
        # the Boresch potential via clambda), so no crgmask.
        if self._stage == "restraint":
            cm = getattr(self.cfg, "restraint_crgmask", "")
            if cm:
                return f"  crgmask='{cm}',\n"
            return ""
        if self._stage == "vdw":
            cm = getattr(self.cfg, "crgmask", ":LIG")
            return f"  crgmask='{cm}',\n"
        return ""

    @property
    def _active_lambdas(self):
        if self._stage == "restraint":
            return list(
            getattr(self.cfg, "restraint_lambdas_fine", None)
            if getattr(self.cfg, "use_fine_restraint_lambdas", False)
            else getattr(self.cfg, "restraint_lambdas", self.cfg.lambdas)
        )
        if self._stage == "decharge":
            return list(getattr(self.cfg, "decharge_lambdas", self.cfg.lambdas))
        if self._stage == "vdw":
            if getattr(self.cfg, "use_dense_vdw", False) and \
                    getattr(self.cfg, "vdw_lambdas_dense", None):
                return list(self.cfg.vdw_lambdas_dense)
            return list(getattr(self.cfg, "vdw_lambdas", self.cfg.lambdas))
        return list(self.cfg.lambdas)

    def _posres_block(self, which):
        """v2.5.19: positional-restraint &cntrl lines for equilibration.

        Holds solute heavy atoms (and, for the complex, the receptor) while the
        box and soft-core ligand relax, so the ligand cannot drift before
        production (FEP-SPell staged-equilibration behaviour; patch sec.4).
        Returns '' when no mask/weight is configured (-> ntr=0 as before).
        `which` selects the weight: 'dens' uses min_restraint_wt (stiff),
        'eq' uses equil_restraint_wt (gentle).
        """
        mask = (getattr(self.cfg, "posres_mask", "") or
                getattr(self.cfg, "posres_mask_default", "") or "")
        if not mask:
            return ""
        if which == "dens":
            wt = getattr(self.cfg, "min_restraint_wt", 5.0)
        else:
            wt = getattr(self.cfg, "equil_restraint_wt", 1.0)
        return (f"  ntr=1, restraint_wt={wt},\n"
                f"  restraintmask='{mask}',\n")

    def _vlimit_block(self):
        v = getattr(self.cfg, "vlimit", 0.0) or 0.0
        return f"  vlimit={v},\n" if v > 0 else ""

    def _vlimit_block_floor(self, floor=20.0):
        """v2.5.59: the restraint-leg prod path was the ONLY production stage
        without a vlimit guard. Four mid-attachment windows (lambda 0.10-0.26)
        suffered single-step velocity explosions (BOND -> 1e8, VOLUME -20%% in
        one 2 ps block; energies perfectly healthy in the prior block ->
        instantaneous clash, NOT a drift/timestep problem). vlimit caps that
        runaway step. Use cfg.vlimit if set higher, else a hard floor."""
        v = getattr(self.cfg, "vlimit", 0.0) or 0.0
        v = max(v, floor)
        return f"  vlimit={v},\n"

    def _ewald_block(self):
        """v2.5.58: skinnb belongs in the &ewald namelist (NOT &cntrl).
        pmemd default skinnb=2.0; raising to 3.0 enlarges the nonbond pairlist
        skin so a benign NPT volume fluctuation does not trigger the spurious
        'Periodic box dimensions have changed too much' GPU cell-list abort.
        Verified against AMBER manual / ParmEd ewald namelist + AMBER devs
        list (fix = add &ewald with skinnb=3.0)."""
        return "&ewald\n  skinnb=3.0,\n /\n"

    def _heat_in(self, cl):
        """v2.5.25/28: NVT heat (ntb=1,ntp=0) with TEMP0 ramp + posres, before the
        first barostat. Single nmropt + single &wt END; DISANG when Boresch on.
        v2.5.40: TI keywords routed through _ti_kw_line so the Option A restraint
        stage (icfe=0, plain MD) emits NONE -- fixes the heat-stage
        'timask1/2 must match the same number of atoms' abort."""
        _dt   = getattr(self.cfg, "heat_dt_ps", 0.001)
        _nst  = getattr(self.cfg, "heat_nstlim", 100_000)
        _gln  = getattr(self.cfg, "heat_gamma_ln", 5.0)
        _t0   = self.cfg.temperature_K
        _ti   = getattr(self.cfg, "heat_T_start", 5.0)
        _wt   = getattr(self.cfg, "heat_restraint_wt", 5.0)
        _mask = (getattr(self.cfg, "posres_mask", "") or
                 getattr(self.cfg, "posres_mask_default", "") or "")
        _ramp = int(_nst * getattr(self.cfg, "heat_ramp_frac", 0.8))
        _pr = (f"  ntr=1, restraint_wt={_wt},\n  restraintmask='{_mask}',\n"
               if _mask else "  ntr=0,\n")
        _disang = "DISANG=boresch.RST\n" if self.boresch else ""
        return f"""FEP NVT heat (TEMP0 ramp), lambda={cl}
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim={_nst}, dt={_dt},
  ntc=2, ntf=1, ntt=3, gamma_ln={_gln}, ig=-1,
  tempi={_ti}, temp0={_t0},
  ntb=1, ntp=0, cut={self.cfg.cutoff_A},
{_pr}  nscm=0, ntwr={_nst // 4},
  ntpr={self.cfg.ntpr}, ntwx=0, ioutfm=1,
{self._vlimit_block()}{self._ti_kw_line(cl)}  nmropt=1,
 /
 &wt
 type = 'TEMP0',
 istep1 = 0, istep2 = {_ramp},
 value1 = {_ti}, value2 = {_t0},
 /
 &wt
 type = 'END',
 /
{_disang}"""

    def _min_in(self, cl):
        rcntrl, rtail = self._restraint_block()
        t1, t2, s1, s2 = self._mask_block()
        return f"""FEP minimization, lambda={cl}
&cntrl
  imin=1, ntmin=2, maxcyc=5000, ncyc=2500,
  ntb=1, ntr=0, cut={self.cfg.cutoff_A},
{self._ti_kw_line(cl)}{rcntrl}/
{rtail}"""

    def _eq_in(self, cl):
        # final60 BUG 1 FIX: equilibration previously jumped straight back to
        # dt=dt_ps (2 fs) + MC barostat off the GENTLE dens.rst, restarting at
        # tempi=100 K with a weak thermostat (gamma_ln=2). On a half-decoupled
        # soft-core ligand at mid-lambda (e.g. complex_vdw lambda=0.50) that is
        # exactly the regime that diverges -- observed T=14,963 K, BOND=3.8M.
        # Fix: inherit the gentle integrator for equilibration too --
        #   * dt = eq_dt_ps (default 0.001, half step)
        #   * Berendsen barostat (eq_barostat=1) with loose taup (eq_taup=5.0)
        #   * stronger Langevin coupling (eq_gamma_ln=5.0)
        #   * restart AT temp0 (tempi=temp0), not 100 K, since dens already
        #     equilibrated temperature.
        # Production keeps the MC barostat / 2 fs step (see _prod_in).
        rcntrl, rtail = self._restraint_block()
        t1, t2, s1, s2 = self._mask_block()
        _eq_dt   = getattr(self.cfg, "eq_dt_ps", 0.001)
        _eq_gln  = getattr(self.cfg, "eq_gamma_ln", 5.0)
        _eq_taup = getattr(self.cfg, "eq_taup", 5.0)
        _eq_bar  = getattr(self.cfg, "eq_barostat", 1)
        # keep the same physical equilibration time despite the smaller step.
        _eq_steps = int(self.cfg.nstlim_eq * (self.cfg.dt_ps / _eq_dt))
        return f"""FEP equilibration NPT, lambda={cl}
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim={_eq_steps}, dt={_eq_dt},
  ntc=2, ntf=1, ntt=3, gamma_ln={_eq_gln}, ig=-1,
  tempi={self.cfg.temperature_K}, temp0={self.cfg.temperature_K},
  ntp=1, pres0=1.0, taup={_eq_taup}, barostat={_eq_bar}, ntb=2, cut={self.cfg.cutoff_A},
{self._vlimit_block()}{self._posres_block('eq')}  nscm=0, ntwr={_eq_steps // 5},
  ntpr={self.cfg.ntpr}, ntwx={self.cfg.ntwx},
{self._ti_kw_line_mbar(cl)}{rcntrl}/
{self._ewald_block()}{rtail}"""

    def _dens_in(self, cl):
        # v2.4.26: GENTLE density pre-equilibration. Old protocol (taup=1.0,
        # dt=0.002, MC barostat from minimized coords) over-corrected the box
        # around a half-decoupled soft-core ligand -> volume doubled, density
        # ~0.44, SC_VDW overflow -> NaN at lambda 0.60/0.65. Fix: taup=5.0
        # (loose), dt=0.001 (half step), Berendsen barostat (barostat=1),
        # stronger thermostat (gamma_ln=5.0). Production keeps MC barostat.
        rcntrl, rtail = self._restraint_block()
        t1, t2, s1, s2 = self._mask_block()
        _dt_dens = 0.001
        _ndens = max(100000, self.cfg.nstlim_eq)
        _taup = getattr(self.cfg, "dens_taup", 5.0)
        _barostat = getattr(self.cfg, "dens_barostat", 1)
        return f"""FEP density pre-equilibration NPT, lambda={cl}
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim={_ndens}, dt={_dt_dens},
  ntc=2, ntf=1, ntt=3, gamma_ln=5.0, ig=-1,
  tempi=100.0, temp0={self.cfg.temperature_K},
  ntp=1, pres0=1.0, taup={_taup}, barostat={_barostat}, ntb=2, cut={self.cfg.cutoff_A},
{self._vlimit_block()}{self._posres_block('dens')}  nscm=0, ntwr={_ndens // 4},
  ntpr={self.cfg.ntpr}, ntwx={self.cfg.ntwx * 10},
{self._ti_kw_line_mbar(cl)}{rcntrl}/
{self._ewald_block()}{rtail}"""

    def _prod_in(self, cl):
        rcntrl, rtail = self._restraint_block()
        if self._stage == "restraint":
            return self._prod_in_restraint(cl, rcntrl, rtail)
        t1, t2, s1, s2 = self._mask_block()
        if self.hremd:
            nchunks = max(1, self.cfg.nstlim_prod // self.exchange_freq)
            rem_block = (f"  numexchg={nchunks},\n"
                         f"  nstlim={self.exchange_freq},\n")
        else:
            rem_block = f"  nstlim={self.cfg.nstlim_prod},\n"
        return f"""FEP production NPT, lambda={cl}
&cntrl
  imin=0, irest=1, ntx=5,
{rem_block}  dt={self.cfg.dt_ps},
  ntc=2, ntf=1, ntt=3, gamma_ln=2.0, ig=-1, temp0={self.cfg.temperature_K},
  ntp=1, pres0=1.0, taup=5.0, barostat=1, ntb=2, cut={self.cfg.cutoff_A},
{self._vlimit_block()}  nscm=0,
  ntpr={self.cfg.ntpr}, ntwx={self.cfg.ntwx}, ioutfm=1,
{self._ti_kw_line_mbar(cl)}{rcntrl}/
{self._ewald_block()}{rtail}"""

    def _prod_in_restraint(self, cl, rcntrl, rtail):
        """v2.5.36 (Option A) restraint-leg production mdin: plain NPT MD.

        v2.5.64: uses cfg.restraint_nstlim_prod (default 2 ns @1fs) instead of
        the full TI nstlim_prod -- this leg is equilibration for an analytic
        Boresch correction, not a sampling-limited TI/MBAR leg. Falls back to
        nstlim_prod if restraint_nstlim_prod is unset/0.
        """
        _rnst = getattr(self.cfg, "restraint_nstlim_prod", 0) or self.cfg.nstlim_prod
        if self.hremd:
            nchunks = max(1, _rnst // self.exchange_freq)
            rem_block = (f"  numexchg={nchunks},\n  nstlim={self.exchange_freq},\n")
        else:
            rem_block = f"  nstlim={_rnst},\n"
        return f"""FEP restraint production NPT, lambda={cl}
&cntrl
  imin=0, irest=1, ntx=5,
{rem_block}  dt={self.cfg.dt_ps},
  ntc=2, ntf=1, ntt=3, gamma_ln=2.0, ig=-1, temp0={self.cfg.temperature_K},
  ntp=1, pres0=1.0, taup=5.0, barostat=1, ntb=2, cut={self.cfg.cutoff_A},
{self._vlimit_block_floor()}  nscm=0,
  ntpr={self.cfg.ntpr}, ntwx={self.cfg.ntwx}, ioutfm=1,
{rcntrl}/
{self._ewald_block()}{rtail}"""

    def _write_boresch_RST(self, wd, lam=1.0):
        """Write the six Boresch &rst records in FEP-SPell-ABFE atom ordering.

        v2.5.18 CORRECTNESS FIX: the previous writer restrained a different
        set of angles/dihedrals than the workflow measured and analytically
        corrected (only 1/6 DOF matched). The simulated restraint now uses
        exactly the FEP-SPell coordinates:
            r=L1-P1, alpha=P1-L1-L2, theta=P2-P1-L1,
            gamma=P1-L1-L2-L3, beta=P2-P1-L1-L2, phi=P3-P2-P1-L1
        with L1=A,L2=B,L3=C (ligand) and P1=aA,P2=bA,P3=cA (receptor).
        """
        b = self.boresch
        # canonical 1-based atom indices
        L1 = int(b.get("L1", b["A"]));  L2 = int(b.get("L2", b["B"]));  L3 = int(b.get("L3", b["C"]))
        P1 = int(b.get("P1", b["aA"])); P2 = int(b.get("P2", b["bA"])); P3 = int(b.get("P3", b["cA"]))
        r0     = float(b["r0"])
        alpha0 = float(b.get("alpha0", b["thB0"]))   # P1-L1-L2
        theta0 = float(b.get("theta0", b["thA0"]))   # P2-P1-L1
        gamma0 = float(b.get("gamma0", b["phC0"]))   # P1-L1-L2-L3
        beta0  = float(b.get("beta0",  b["phB0"]))   # P2-P1-L1-L2
        phi0   = float(b.get("phi0",   b["phA0"]))   # P3-P2-P1-L1
        # v2.5.32: restraint-leg = Boresch force-constant SCALING.
        #   * restraint stage: rk(lambda) = lambda * k  (0 at lambda=0 -> full at
        #     lambda=1). Both TI ends are the same fully-interacting complex, so
        #     the ONLY lambda dependence is the restraint -> small, smooth dV/dl.
        #   * decharge / vdw stages: FIXED full k (restraint HELD ON the whole
        #     time so the ligand stays in the pocket while it is decoupled).
        # v2.5.33: FIXED Boresch force constants for ALL stages/legs (FEP-SPell
        # recipe). The restraint-leg lambda dependence comes from the :1<->:2
        # dual-copy TI transformation, NOT from scaling the restraint.
        _scale = float(lam)
        kr  = float(b["kr"])  * _scale
        kth = float(b["kth"]) * _scale
        kph = float(b["kph"]) * _scale
        lines = [_RESTRAINT_HEADER]
        # distance r : L1-P1
        lines.append(self._rst_line([L1, P1],
            0.0, r0, r0, 99.0, kr, kr))
        # angle alpha : P1-L1-L2
        lines.append(self._rst_line([P1, L1, L2],
            0.0, alpha0, alpha0, 180.0, kth, kth))
        # angle theta : P2-P1-L1
        lines.append(self._rst_line([P2, P1, L1],
            0.0, theta0, theta0, 180.0, kth, kth))
        # dihedral gamma : P1-L1-L2-L3
        lines.append(self._rst_line([P1, L1, L2, L3],
            gamma0-180.0, gamma0, gamma0, gamma0+180.0, kph, kph))
        # dihedral beta : P2-P1-L1-L2
        lines.append(self._rst_line([P2, P1, L1, L2],
            beta0-180.0, beta0, beta0, beta0+180.0, kph, kph))
        # dihedral phi : P3-P2-P1-L1
        lines.append(self._rst_line([P3, P2, P1, L1],
            phi0-180.0, phi0, phi0, phi0+180.0, kph, kph))
        # v2.5.31d: refuse to emit a restraint whose reference does not match coords
        self._verify_boresch_reference(wd, L1, P1, r0)
        (wd / "boresch.RST").write_text("\n".join(lines) + "\n")


    @staticmethod
    def _read_rst7_coords_box(rst7):
        """Parse an ASCII Amber rst7: return (coords list[(x,y,z)], box or None)."""
        try:
            raw = open(rst7).read().splitlines()
        except Exception:
            return None, None
        if len(raw) < 2:
            return None, None
        try:
            natom = int(raw[1].split()[0])
        except Exception:
            return None, None
        nums = []
        ncoord_lines = (natom * 3 + 5) // 6
        for ln in raw[2:2 + ncoord_lines]:
            for j in range(0, len(ln), 12):
                tok = ln[j:j + 12].strip()
                if tok:
                    try:
                        nums.append(float(tok))
                    except ValueError:
                        pass
        if len(nums) < natom * 3:
            return None, None
        coords = [(nums[3 * k], nums[3 * k + 1], nums[3 * k + 2]) for k in range(natom)]
        box = None
        last = raw[-1].split()
        if len(last) >= 6:
            try:
                box = tuple(float(x) for x in last[:6])
            except ValueError:
                box = None
        return coords, box

    @staticmethod
    def _min_image_dist(a, b, box):
        """L1-P1 distance with orthorhombic minimum image (oct approximated by lengths)."""
        import math
        dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
        if box:
            lx, ly, lz = box[0], box[1], box[2]
            if lx > 0:
                dx -= lx * round(dx / lx)
            if ly > 0:
                dy -= ly * round(dy / ly)
            if lz > 0:
                dz -= lz * round(dz / lz)
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def _verify_boresch_reference(self, wd, L1, P1, r0):
        """v2.5.31d ROOT-CAUSE GUARD: fail loudly at WRITE time if the Boresch
        reference r0 does not match the coordinates the leg starts from, or if
        the anchors span a periodic image. Such a mismatch injects ~k*delta^2
        kcal/mol and aborts GPU production with 'Periodic box dimensions changed
        too much' after wasted restarts."""
        for cand in ("system.inpcrd", "min.rst", "dens.rst"):
            crd = wd / cand
            if not crd.exists():
                continue
            coords, box = self._read_rst7_coords_box(crd)
            if not coords or max(L1, P1) > len(coords):
                continue
            actual = self._min_image_dist(coords[L1 - 1], coords[P1 - 1], box)
            half = min(box[0], box[1], box[2]) / 2.0 if box else None
            problems = []
            if abs(actual - r0) > 2.0:
                problems.append(
                    "actual L1-P1 distance %.2f A vs reference r0=%.2f A (delta %+.2f)"
                    % (actual, r0, actual - r0))
            if half and actual > half:
                problems.append(
                    "L1-P1 distance %.2f A exceeds half-box %.2f A (periodic-image anchor)"
                    % (actual, half))
            if problems:
                raise ValueError(
                    "Boresch reference geometry does not match start coordinates in "
                    + str(wd) + ": " + "; ".join(problems)
                    + ". Re-derive Boresch anchors/r0 from the SAME frame the leg "
                      "starts from (and apply minimum image).")
            return
        return

    @staticmethod
    def _rst_line(iat, r1, r2, r3, r4, rk2, rk3):
        atoms = ",".join(str(i) for i in iat) + ","
        return (f"&rst iat={atoms} "
                f"r1={r1:.3f}, r2={r2:.3f}, r3={r3:.3f}, r4={r4:.3f}, "
                f"rk2={rk2:.3f}, rk3={rk3:.3f}, /")

    def setup_leg(self, leg, prmtop, inpcrd, stage=None,
                  write_correction=True):
        self._stage = stage
        self._suppress_correction = not write_correction
        leg_dir = ensure_dir(self.root / leg)
        for lam in self._active_lambdas:
            wd = ensure_dir(leg_dir / f"lambda_{lam:.3f}")
            (wd / "min.in").write_text(self._min_in(lam))
            if getattr(self.cfg, "do_heat", True):
                (wd / "heat.in").write_text(self._heat_in(lam))
            (wd / "dens.in").write_text(self._dens_in(lam))   # v2.4.22
            (wd / "eq.in").write_text(self._eq_in(lam))
            (wd / "prod.in").write_text(self._prod_in(lam))
            if self.boresch:
                self._write_boresch_RST(wd, lam=lam)
                _bm = getattr(self, "boresch", None) or {}
                _lm = _bm.get("lig_masks"); _rm = _bm.get("rec_masks")
                if _lm and _rm:
                    import json as _json
                    (wd / "boresch_eqcheck.json").write_text(
                        _json.dumps({"lig_masks": _lm, "rec_masks": _rm}))
            for src, dst in [(prmtop, wd / "system.prmtop"),
                             (inpcrd, wd / "system.inpcrd")]:
                try:
                    if dst.exists() or dst.is_symlink(): dst.unlink()
                    dst.symlink_to(Path(src).resolve())
                except OSError:
                    dst.write_bytes(Path(src).read_bytes())
        if (self.boresch and "dG_correction_kcal_mol" in self.boresch
                and not self._suppress_correction):
            (leg_dir / "boresch_correction.txt").write_text(
                f"{self.boresch['dG_correction_kcal_mol']:.6f}\n")
        # v2.5.42: pre-submission guard -- the Option A restraint leg must be
        # plain MD (icfe=0, no TI keywords) in EVERY stage. If any mdin leaks a
        # TI keyword we fail loudly NOW rather than after a wasted GPU run.
        if self._stage == "restraint":
            self._validate_restraint_inputs(leg_dir)
        return leg_dir

    # v2.5.42: TI keywords that must NEVER appear in an Option A restraint-leg
    # mdin (the leg is single-copy plain MD with a lambda-scaled Boresch &rst).
    _RESTRAINT_FORBIDDEN_KW = ("icfe=1", "timask1", "timask2", "scmask1",
                               "scmask2", "clambda", "ifmbar", "mbar_lambda",
                               "mbar_states", "crgmask")

    def _validate_restraint_inputs(self, leg_dir):
        """Scan every generated mdin in the restraint leg for forbidden TI
        keywords. Raises RuntimeError listing each offending file/keyword so a
        broken Option A conversion (e.g. a missed stage template) is caught
        before submission. Pure read-only check."""
        from pathlib import Path as _P
        offenders = []
        for mdin in sorted(_P(leg_dir).glob("lambda_*/*.in")):
            try:
                text = mdin.read_text()
            except OSError:
                continue
            hits = [kw for kw in self._RESTRAINT_FORBIDDEN_KW if kw in text]
            if hits:
                offenders.append((mdin, hits))
        if offenders:
            lines = ["RESTRAINT-LEG INPUT VALIDATION FAILED (Option A): the "
                     "restraint leg must contain NO TI keywords, but found:"]
            for mdin, hits in offenders:
                rel = mdin.relative_to(leg_dir)
                lines.append(f"  {rel}: {', '.join(hits)}")
            lines.append("This indicates a stage template still emits TI "
                         "keywords. Refusing to stage a run that would abort "
                         "with 'timask1/2 must match' or skip production.")
            raise RuntimeError("\n".join(lines))
        log.info("  [restraint-input-guard] %s: all mdin clean (no TI keywords)",
                 leg_dir.name)
    def _avoid_select_line(self):
        """v2.5.8: build an LSF -R select line excluding hpc_cfg.avoid_hosts so
        every generated GPU job skips known-bad nodes. '' when none set."""
        hosts = sorted(h for h in getattr(self.hpc_cfg, "avoid_hosts", ()) if h)
        if not hosts:
            return ""
        expr = " && ".join(f"hname!='{h}'" for h in hosts)
        return f'#BSUB -R "select[{expr}]"\n'

    def build_lsf_array(self, leg_dir, leg_name):
        nwin = len(self._active_lambdas)
        script = leg_dir / f"run_{leg_name}.lsf"
        h  = "#!/bin/bash\n"
        h += lsf_banner(f"{leg_name} array")
        h += f"#BSUB -q {self.hpc_cfg.queue_gpu}\n"
        h += f"#BSUB -P {self.hpc_cfg.project}\n"
        h += f"#BSUB -J fep_{leg_name}[1-{nwin}]\n"
        h += f"#BSUB -W {self.hpc_cfg.walltime}\n"
        h += f"#BSUB -o fep_{leg_name}.%I.%J.out\n"
        h += f"#BSUB -e fep_{leg_name}.%I.%J.err\n"
        _mem = getattr(self.hpc_cfg, "fep_mem_mb", 8192)
        # v2.5.34: reserve CPU cores alongside the 1 GPU so the rare inline
        # CPU density-settle fallback has parallelism. This builder requests
        # the GPU explicitly via the -gpu "num=" line below, so -n here is
        # purely CPU cores (no GPU-count inflation). GPU MD uses ~1-2 cores;
        # 16 is trivial on 32-64 core the cluster GPU nodes. span[hosts=1] keeps
        # them co-located for the MPI settle. Tune via HPCConfig.fep_gpu_cores.
        _gpu_cores = getattr(self.hpc_cfg, "fep_gpu_cores", 16)
        h += f"#BSUB -n {_gpu_cores}\n"
        h += '#BSUB -R "span[hosts=1]"\n'
        h += f'#BSUB -R "rusage[mem={_mem}]"\n'
        h += self._avoid_select_line()
        h += f'#BSUB -gpu "num={self.hpc_cfg.n_gpu}"\n\n'
        h += "module purge\n"
        for m in self.hpc_cfg.modules: h += f"module load {m}\n"
        if self.hpc_cfg.venv_activate: h += f"source {self.hpc_cfg.venv_activate}\n"
        h += "set -euo pipefail\n\n"
        lam_array = " ".join(f"{l:.3f}" for l in self._active_lambdas)
        # v2.5.31d: values for the GPU regrid step heredoc (f-string scope)
        _rg_t1, _rg_t2, _rg_s1, _rg_s2 = self._mask_block()
        _rg_ifsc = self._ifsc_value()
        _rg_nmropt = 1 if self.boresch else 0
        _rg_tail = "&wt type='END' /\nDISANG=boresch.RST\n" if self.boresch else ""
        body = f"""LAMBDAS=({lam_array})
LAM=${{LAMBDAS[$((LSB_JOBINDEX-1))]}}
WD={leg_dir}/lambda_${{LAM}}
cd "$WD"

# --- v2.4.15: surface the REAL pmemd failure instead of a bare exit 1 ---
run_stage() {{
    # $1 = stage label, remaining args = pmemd command
    local stage="$1"; shift
    local mdout=""
    # the token after -o is this stage's mdout file
    local prev=""
    for a in "$@"; do
        if [ "$prev" = "-o" ]; then mdout="$a"; fi
        prev="$a"
    done
    echo "[run_stage] lambda=${{LAM}} stage=${{stage}} starting: $*"
    set +e
    "$@"
    local rc=$?
    set -e
    if [ "$rc" -ne 0 ]; then
        echo "============================================================" >&2
        echo "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=${{stage}} rc=${{rc}}" >&2
        if [ -n "$mdout" ] && [ -f "$mdout" ]; then
            echo "---- cause (from ${{mdout}}) ----" >&2
            grep -iE 'STOP PMEMD|Terminated Abnormally|NaN|vlimit exceeded|Input errors|ERROR' \\
                 "$mdout" | tail -n 20 >&2 || true
            echo "---- last 15 lines of ${{mdout}} ----" >&2
            tail -n 15 "$mdout" >&2 || true
        else
            echo "(no mdout '$mdout' produced -- pmemd died before writing output)" >&2
        fi
        echo "============================================================" >&2
        exit "$rc"
    fi
}}

run_stage min  pmemd.cuda -O -i min.in  -p system.prmtop -c system.inpcrd \\
               -o min.out  -r min.rst  -x min.nc
# --- v2.5.31k: PRE-FLIGHT Boresch geometry gate (FAIL FAST) -------------------
# Validate the six Boresch atoms against the MINIMISED structure (min.rst),
# BEFORE heat/dens/eq. A broken restraint (e.g. ligand built outside the pocket,
# r~63 A, RESTRAINT~50,000) is caught in seconds here instead of ~50 minutes
# later at the post-eq gate. Same validator; only the structure differs.
#   PASS -> exit 0 -> "preflight geometry OK", continue to heat.
#   FAIL -> exit 73 -> window stops now (no heat/dens/eq wasted).
#   SKIP -> exit 3  -> deps/manifest missing; do NOT fail (post-eq gate still runs).
if [ -f boresch.RST ] && [ -f min.rst ] && command -v python >/dev/null 2>&1; then
    set +e
    python - "$WD" <<'PYPRE'
import sys, json, os
from pathlib import Path
wd = Path(sys.argv[1])
cands = []
if os.environ.get("AMBER_MD_HOME"):
    cands.append(os.environ["AMBER_MD_HOME"]); cands.append(str(Path(os.environ["AMBER_MD_HOME"]).parent))
for p in os.environ.get("PYTHONPATH", "").split(os.pathsep):
    if p: cands.append(p)
here = wd.resolve()
for up in [here] + list(here.parents):
    if (up / "amber_md" / "__init__.py").is_file():
        cands.append(str(up)); break
for c in cands:
    if c and c not in sys.path: sys.path.insert(0, c)
try:
    import parmed as pmd  # noqa: F401
    from amber_md.boresch_autoselect import validate_masks
except Exception as e:
    print("[boresch-preflight] deps unavailable (%s); SKIPPED" % e); sys.exit(3)
bj = wd / "boresch_eqcheck.json"
if not bj.exists():
    print("[boresch-preflight] no eq-check manifest; SKIPPED"); sys.exit(3)
m = json.loads(bj.read_text())
chk = validate_masks(str(wd / "system.prmtop"), str(wd / "min.rst"),
                     m["lig_masks"], m["rec_masks"])
print("[boresch-preflight] r=%.2f alpha=%.1f theta=%.1f -> %s"
      % (chk["r"], chk["alpha"], chk["theta"],
         "PASS" if chk["all_pass"] else "FAIL %s" % chk["failures"]))
sys.exit(0 if chk["all_pass"] else 73)
PYPRE
    pre_rc=$?
    set -e
    if [ "$pre_rc" -eq 0 ]; then
        echo "[boresch-preflight] lambda=${{LAM}} minimised geometry OK"
    elif [ "$pre_rc" -eq 73 ]; then
        echo "============================================================" >&2
        echo "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=preflight: Boresch geometry out of bounds at min.rst (broken restraint -- e.g. ligand not in pocket). Failing now instead of wasting heat/dens/eq. Check boresch.RST iat indices and that build_restraint placed the ligand in the pocket." >&2
        echo "============================================================" >&2
        exit 73
    else
        echo "[boresch-preflight] lambda=${{LAM}} SKIPPED (rc=${{pre_rc}}); post-eq gate will still validate" >&2
    fi
fi
# v2.5.25: explicit NVT heat (ntp=0) before the first barostat
if [ -f heat.in ]; then
run_stage heat pmemd.cuda -O -i heat.in -p system.prmtop -c min.rst -ref min.rst \\
               -o heat.out -r heat.rst -x heat.nc
HEAT_TEMP_MAX_K=${{HEAT_TEMP_MAX_K:-1000}}
if [ -f heat.out ]; then
    maxTh=$(set +o pipefail; set +e; \\
           grep -aoE 'TEMP\\(K\\) *= *[-0-9.]+' heat.out 2>/dev/null \\
           | grep -aoE '[-0-9.]+$' \\
           | awk 'BEGIN{{m=""}} {{if(m==""||$1>m)m=$1}} END{{if(m!="")print m}}'; \\
           true)
    if [ -n "$maxTh" ] && awk -v t="$maxTh" -v m="$HEAT_TEMP_MAX_K" 'BEGIN{{exit !(t>m)}}'; then
        echo "============================================================" >&2
        echo "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=heat INSTABILITY: peak TEMP=${{maxTh}} K > ${{HEAT_TEMP_MAX_K}} K" >&2
        grep -aiE 'vlimit exceeded|NaN|Terminated Abnormally' heat.out | tail -n 20 >&2 || true
        echo "============================================================" >&2
        exit 70
    fi
fi
DENS_C=heat.rst
else
DENS_C=min.rst
fi
run_stage dens pmemd.cuda -O -i dens.in -p system.prmtop -c "$DENS_C" -ref "$DENS_C" \\
               -o dens.out -r dens.rst -x dens.nc
run_stage eq   pmemd.cuda -O -i eq.in   -p system.prmtop -c dens.rst -ref dens.rst \\
               -o eq.out   -r eq.rst   -x eq.nc

# --- final60 BUG 1 FIX: equilibration stability gate ---
# run_stage only checks the exit code; an eq that integrates to 15,000 K but
# still exits 0 would silently pass and poison every downstream estimator.
# Fail the window if the peak temperature in eq.out exceeds EQ_TEMP_MAX_K.
EQ_TEMP_MAX_K=${{EQ_TEMP_MAX_K:-1000}}
if [ -f eq.out ]; then
    # v2.5.7 FIX: compute peak eq temperature WITHOUT a pipeline that 'head'
    # closes early. Under 'set -o pipefail', 'grep ... | sort | head -n1' makes
    # upstream grep/sort die with SIGPIPE (rc=141) once head exits, and
    # 'set -e' then aborts the whole window with exit 141 BEFORE prod runs.
    # That is exactly the all-windows / every-node exit-141 failure. Run the
    # scan in a pipefail/errexit-free subshell; awk does the max in one pass.
    maxT=$(set +o pipefail; set +e; \
           grep -aoE 'TEMP\\(K\\) *= *[-0-9.]+' eq.out 2>/dev/null \
           | grep -aoE '[-0-9.]+$' \
           | awk 'BEGIN{{m=""}} {{if(m==""||$1>m)m=$1}} END{{if(m!="")print m}}'; \
           true)
    if [ -n "$maxT" ]; then
        if awk -v t="$maxT" -v m="$EQ_TEMP_MAX_K" 'BEGIN{{exit !(t>m)}}'; then
            echo "============================================================" >&2
            echo "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=eq INSTABILITY: peak TEMP=${{maxT}} K > ${{EQ_TEMP_MAX_K}} K" >&2
            grep -aiE 'vlimit exceeded|NaN|Terminated Abnormally' eq.out \\
                 | tail -n 20 >&2 || true
            echo "(eq integrated to a non-physical temperature; refusing to " >&2
            echo " run production on a blown-up box. Lower dt/strengthen the " >&2
            echo " thermostat or add lambda windows around this region.)" >&2
            echo "============================================================" >&2
            exit 71
        fi
    fi
fi
# --- v2.5.17 (B): post-equilibration Boresch geometry GATE --------------------
# Re-validate the SAME six Boresch atoms against the EQUILIBRATED structure
# (eq.rst), AFTER the eq temperature gate and BEFORE production. Outcomes:
#   PASS  -> exit 0  -> prints "...geometry OK", production proceeds.
#   FAIL  -> exit 72 -> window stops (geometry drifted out of bounds).
#   SKIP  -> exit 3  -> prints "...SKIPPED (NOT validated)"; production still
#            proceeds, but NO false "OK" is printed (v2.5.17 import-robust fix).
# Import robustness: the heredoc runs a bare interpreter whose sys.path may not
# contain the amber_md install. We inject candidate roots (AMBER_MD_HOME,
# PYTHONPATH entries, and parents of $WD) before importing.
if [ -f boresch.RST ] && [ -f eq.rst ] && command -v python >/dev/null 2>&1; then
    set +e
    python - "$WD" <<'PYGATE'
import sys, json, os
from pathlib import Path
wd = Path(sys.argv[1])

# --- make `import amber_md` work on the compute node ---------------------------
cands = []
if os.environ.get("AMBER_MD_HOME"):
    cands.append(os.environ["AMBER_MD_HOME"])
    cands.append(str(Path(os.environ["AMBER_MD_HOME"]).parent))
for p in os.environ.get("PYTHONPATH", "").split(os.pathsep):
    if p:
        cands.append(p)
# walk up from the lambda working dir: .../<leg>/lambda_x/ -> repo root that
# contains an `amber_md/` package directory.
here = wd.resolve()
for up in [here] + list(here.parents):
    if (up / "amber_md" / "__init__.py").is_file():
        cands.append(str(up)); break
for c in cands:
    if c and c not in sys.path:
        sys.path.insert(0, c)

try:
    import parmed as pmd  # noqa: F401
    from amber_md.boresch_autoselect import validate_masks
except Exception as e:
    print("[boresch-gate] deps unavailable (%s); SKIPPED (NOT validated)" % e)
    sys.exit(3)

bj = wd / "boresch_eqcheck.json"
if not bj.exists():
    print("[boresch-gate] no eq-check manifest; SKIPPED (NOT validated)")
    sys.exit(3)

m = json.loads(bj.read_text())
chk = validate_masks(str(wd / "system.prmtop"), str(wd / "eq.rst"),
                     m["lig_masks"], m["rec_masks"])
print("[boresch-gate] r=%.2f alpha=%.1f theta=%.1f -> %s"
      % (chk["r"], chk["alpha"], chk["theta"],
         "PASS" if chk["all_pass"] else "FAIL %s" % chk["failures"]))
sys.exit(0 if chk["all_pass"] else 72)
PYGATE
    gate_rc=$?
    set -e
    if [ "$gate_rc" -eq 0 ]; then
        echo "[boresch-gate] lambda=${{LAM}} equilibrated geometry OK"
    elif [ "$gate_rc" -eq 72 ]; then
        echo "[boresch-gate] lambda=${{LAM}} FAILED: equilibrated Boresch geometry out of bounds; refusing production" >&2
        exit 72
    elif [ "$gate_rc" -eq 3 ]; then
        echo "[boresch-gate] lambda=${{LAM}} FAILED: validation deps/manifest unavailable on node (amber_md not importable). Refusing to run production UNVALIDATED. Install the package on the compute node (pip install -e .)." >&2
        exit 3
    else
        echo "[boresch-gate] lambda=${{LAM}} FAILED: gate error rc=${{gate_rc}}; refusing production" >&2
        exit "$gate_rc"
    fi
fi
# --- v2.4.22: auto-restart prod on GPU "box changed too much" (rc=255) ---
# pmemd.cuda halts (not NaN) when NPT volume drift invalidates the fixed
# GPU grid. Amber's prescribed fix is to restart from the latest restart
# file, which rebuilds the grid. We retry up to 5 times, each time feeding
# the most recent restart (prod.rst once it exists, else good_rst). NOTE:
# good_rst is initialised from eq.rst but is OVERWRITTEN by each successful
# regrid, so the variable name "eq.rst" is the file path, not necessarily the
# original equilibrated coords -- orig_eq.rst (v2.5.49) preserves those.
# prod.in uses irest=1/ntx=5 so every restart continues seamlessly.
run_prod_with_restart() {{
    local max=10 attempt=1 rc=0
    local good_rst="eq.rst"
    local cpu_used=0
    # v2.5.49: keep a PRISTINE copy of the true equilibrated restart. The regrid
    # step overwrites good_rst (named eq.rst) with settled coords, so without this
    # we permanently lose the real eq geometry. orig_eq.rst lets an escalated
    # settle restart from true equilibrium instead of re-seeding the same clash.
    if [ ! -s orig_eq.rst ] && [ -s eq.rst ]; then cp -f eq.rst orig_eq.rst; fi
    local _drift_count=0
    local _clash_rewound=0   # v2.5.60: ensures the origin-rewind happens at most once
    # v2.5.35: durable recovery log (LSF .err is unreliable on the cluster array tasks)
    local REC_LOG="recovery.log"
    _rec() {{ echo "$@" | tee -a "$REC_LOG" >&2; }}
    # Parse immutable physics from prod.in ONCE so recovery inputs are built
    # from known-good values rather than by mutating the file in place.
    # v2.5.44: each parse ends with '|| true'. Under 'set -euo pipefail' a
    # grep that finds NO match (e.g. Option A prod.in has no clambda/timask/
    # ifsc/icfe) exits 1; a BARE 'VAR=$(failing-pipeline)' assignment would
    # then abort the whole window BEFORE production -- the silent 'eq done but
    # no prod.out' failure. '|| true' makes a no-match yield an empty string.
    local _CL _T0 _TM1 _TM2 _IFSC _NMR _ICFE
    _CL=$(grep -oE 'clambda=[0-9.]+' prod.in | head -1 | cut -d= -f2) || true
    _T0=$(grep -oE 'temp0=[0-9.]+'  prod.in | head -1 | cut -d= -f2) || true
    _TM1=$(grep -oE "timask1='[^']*'" prod.in | head -1 | cut -d= -f2) || true
    _TM2=$(grep -oE "timask2='[^']*'" prod.in | head -1 | cut -d= -f2) || true
    _IFSC=$(grep -oE 'ifsc=[0-9]+' prod.in | head -1 | cut -d= -f2) || true
    _ICFE=$(grep -oE 'icfe=[0-9]+' prod.in | head -1 | cut -d= -f2) || true
    _NMR=$(grep -oE 'nmropt=[0-9]+' prod.in | head -1 | cut -d= -f2) || true
    : "${{_T0:=298.0}}"; : "${{_IFSC:=0}}"; : "${{_NMR:=0}}"; : "${{_ICFE:=0}}"
    # v2.5.41: the Option A restraint leg is plain MD (icfe=0) and has NO
    # clambda -- correct, not an error. Only a TI leg (icfe=1) needs clambda.
    if [ "$_ICFE" = "1" ] && [ -z "$_CL" ]; then
        _rec "[run_stage] lambda=${{LAM}} ABORT: icfe=1 but no clambda in prod.in"
        return 74
    fi
    : "${{_CL:=0.0}}"
    # $1=outfile $2=nstlim $3=dt $4=taup $5=label
    _write_recovery_mdin() {{
        local _out="$1" _nst="$2" _dt="$3" _taup="$4" _lab="$5"
        {{
            echo "FEP $_lab NPT, lambda=$_CL"
            echo "&cntrl"
            echo "  imin=0, irest=0, ntx=1,"
            echo "  nstlim=$_nst, dt=$_dt,"
            echo "  ntc=2, ntf=1, ntt=3, gamma_ln=2.0, ig=-1,"
            echo "  tempi=$_T0, temp0=$_T0,"
            echo "  ntp=1, barostat=1, pres0=1.0, taup=$_taup, ntb=2, cut=10.0,"
            echo "  vlimit=20.0, nscm=0,"
            echo "  ntpr=1000, ntwx=0, ioutfm=1,"
            # v2.5.41: only a TI leg (icfe=1) gets icfe/ifsc/clambda + timask.
            if [ "$_ICFE" = "1" ]; then
            echo "/"
            echo "&ewald"
            echo "  skinnb=3.0,"
            echo " /"
                if [ -n "$_TM1" ]; then echo "  timask1=$_TM1, timask2=$_TM2,"; fi
            fi
            if [ "$_NMR" = "1" ]; then echo "  nmropt=1,"; fi
            echo "/"
            if [ "$_NMR" = "1" ]; then
                echo "&wt type='END' /"
                echo "DISANG=boresch.RST"
            fi
        }} > "$_out"
    }}
    local _alloc="${{LSB_DJOB_NUMPROC:-1}}"
    local _cap=48
    if [ -n "${{CUDA_VISIBLE_DEVICES:-}}" ] || [ -n "${{LSB_GPU_ALLOC:-}}" ] \
       || echo "${{LSB_QUEUE:-}}" | grep -qiE 'gpu'; then
        _cap=32
    fi
    local CPU_NP="${{AMBERMD_CPU_NP:-$_alloc}}"
    if [ "$CPU_NP" -gt "$_cap" ]; then CPU_NP="$_cap"; fi
    if [ "$CPU_NP" -lt 1 ]; then CPU_NP=1; fi
    local CPU_CMD=""
    if command -v pmemd.MPI >/dev/null 2>&1 && command -v mpirun >/dev/null 2>&1 && [ "$CPU_NP" -gt 1 ]; then
        CPU_CMD="mpirun --mca btl_tcp_if_include lo --mca orte_base_help_aggregate 0 -np ${{CPU_NP}} pmemd.MPI"
    elif command -v pmemd >/dev/null 2>&1; then
        CPU_CMD="pmemd"
    elif command -v sander >/dev/null 2>&1; then
        CPU_CMD="sander"
    fi
    cpu_density_settle() {{
        if [ -z "$CPU_CMD" ]; then
            _rec "[CPU_FALLBACK] lambda=${{LAM}} UNAVAILABLE (no pmemd.MPI/pmemd/sander on PATH)"
            return 1
        fi
        _rec "============================================================"
        _rec "[CPU_FALLBACK] lambda=${{LAM}} ENGAGED: GPU box-drift recovery exhausted;"
        _rec "[CPU_FALLBACK] lambda=${{LAM}} running ~25 ps CPU density-settle via: ${{CPU_CMD}}"
        _rec "[CPU_FALLBACK] lambda=${{LAM}} (alloc=${{_alloc}}, cap=${{_cap}}, using -np=${{CPU_NP}}, queue=${{LSB_QUEUE:-?}}); host=$(hostname)"
        echo "============================================================" >&2
        # v2.5.35: clean-template settle input (no fragile sed on prod.in)
        _write_recovery_mdin prod_cpusettle.in 25000 0.001 2.0 "CPU density-settle"
        set +e
        ${{CPU_CMD}} -O -i prod_cpusettle.in -p system.prmtop -c "$good_rst" \
                     -o prod_cpusettle.out -r prod_cpusettle.rst -x prod_cpusettle.nc \
                     > prod_cpusettle.console 2>&1
        local crc=$?; set -e
        if [ "$crc" -eq 0 ] && [ -s prod_cpusettle.rst ]; then
            cp -f prod_cpusettle.rst "$good_rst"; rm -f prod.rst; cpu_used=1
            _rec "[CPU_FALLBACK] lambda=${{LAM}} SUCCESS: box density settled on CPU; resuming GPU prod"
            return 0
        fi
        _rec "[CPU_FALLBACK] lambda=${{LAM}} FAILED (rc=${{crc}}); see prod_cpusettle.console"
        tail -n 25 prod_cpusettle.console >&2 || true; return 1
    }}
    while [ "$attempt" -le "$max" ]; do
        local rsrc="$good_rst"
        if [ -s prod.rst ]; then rsrc="prod.rst"; fi
        echo "[run_stage] lambda=${{LAM}} prod attempt ${{attempt}}/${{max}} from ${{rsrc}}"
        set +e
        pmemd.cuda -O -i prod.in -p system.prmtop -c "$rsrc" \
                   -o prod.out -r prod.rst -x prod.nc > prod.console.${{attempt}} 2>&1
        rc=$?; set -e
        cat prod.console.${{attempt}} || true
        if [ "$rc" -eq 0 ]; then
            # v2.5.59: NaN/explosion sanity gate. A window can exit rc=0 yet have
            # detonated (single-step velocity blow-up: TEMP=NaN, BOND~1e8). Such a
            # prod.out must NOT be reported OK -- it would silently poison MBAR/TI.
            # v2.5.70 FIX: scope the gate to the MD TRAJECTORY ONLY. The
            # 'MBAR Energy analysis:' tables legitimately contain '****' field
            # overflows (Fortran width) and huge finite cross-terms at DISTANT
            # off-diagonal lambda pairs -- these are handled downstream by
            # _sanitize_u_nk and are NOT a detonation. Grepping the whole file
            # false-flagged HEALTHY high-lambda windows (0.725-1.000) as
            # 'rc=0-but-NaN' and discarded good production. We therefore strip
            # everything from the first 'MBAR Energy analysis:' line onward and
            # check only the real per-step MD energy records.
            sed '/MBAR Energy analysis:/,$d' prod.out > .prod_mdonly.$$ 2>/dev/null || cp prod.out .prod_mdonly.$$
            if grep -qiE 'TEMP\(K\) =[[:space:]]*NaN|Etot[[:space:]]*=[[:space:]]*NaN|\*\*\*\*\*\*\*\*\*\*' .prod_mdonly.$$ 2>/dev/null; then
                _rec "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=prod rc=0-but-NaN (MD trajectory detonated despite vlimit; refusing to report OK)"
                grep -niE 'TEMP\(K\) =[[:space:]]*NaN|Etot[[:space:]]*=[[:space:]]*NaN' .prod_mdonly.$$ 2>/dev/null | head -n 3 >&2 || true
                rm -f .prod_mdonly.$$ || true
                return 255
            fi
            rm -f .prod_mdonly.$$ || true
            _rec "[run_stage] lambda=${{LAM}} prod completed OK (attempt ${{attempt}})"; return 0
        fi
        if grep -qiE 'box dimensions have changed too much|changed too much from their initial|Periodic box dimensions have changed' \
                prod.console.${{attempt}} prod.out 2>/dev/null; then
            # v2.5.52 FIX2: accept a drift-halt that occurred AT/NEAR completion.
            # pmemd runs nstlim then can drift-halt on the final box check; the
            # sampling is already done, so do NOT rewind the whole segment.
            _tgt_nst=$(grep -oE 'nstlim=[0-9]+' prod.in | head -1 | cut -d= -f2) || true
            _last_nst=$(grep -oE 'NSTEP =[[:space:]]*[0-9]+' prod.out 2>/dev/null | tail -1 | grep -oE '[0-9]+$') || true
            if [ -n "${{_tgt_nst:-}}" ] && [ -n "${{_last_nst:-}}" ] \
               && [ "$_last_nst" -ge "$(( _tgt_nst * 9 / 10 ))" ]; then
                # v2.5.53: step-count is necessary but NOT sufficient. Require the
                # run to also be PHYSICALLY HEALTHY at the end before accepting it
                # (don't replace a trigger-happy guard with a blind one). Checks:
                #   (a) no NaN anywhere in prod.out
                #   (b) last reported Density in a physical range [0.85,1.15] g/cm^3
                _has_nan=0
                if grep -qiE 'NaN' prod.out 2>/dev/null; then _has_nan=1; fi
                _last_dens=$(grep -oE 'Density[[:space:]]*=[[:space:]]*[0-9.]+' prod.out 2>/dev/null | tail -1 | grep -oE '[0-9.]+$') || true
                _dens_ok=0
                if [ -n "${{_last_dens:-}}" ]; then
                    _dens_ok=$(awk -v d="$_last_dens" 'BEGIN{{print (d>=0.85 && d<=1.15)?1:0}}')
                fi
                if [ "$_has_nan" -eq 0 ] && [ "$_dens_ok" -eq 1 ]; then
                    _rec "[run_stage] lambda=${{LAM}} box-drift halt at completion (NSTEP=${{_last_nst}}/${{_tgt_nst}} >=90%, density=${{_last_dens}}, no NaN); accepting prod as converged (no rewind)"
                    if [ -s prod.rst ]; then cp -f prod.rst "$good_rst"; fi
                    return 0
                else
                    _rec "[run_stage] lambda=${{LAM}} box-drift halt at >=90% but FAILED health gate (nan=${{_has_nan}}, density=${{_last_dens:-none}}); NOT accepting -- continuing recovery ladder"
                fi
            fi
            _drift_count=$((_drift_count+1))
            # v2.5.49: ESCALATING regrid ladder. A short 50 ps @ 1 fs regrid clears
            # slow NPT volume drift (most windows), but a persistent STERIC CLASH
            # (e.g. lambda=0.150: VDW spike + PRESS>8000 at the same step every
            # time) is not resolved by repeating the same gentle regrid from the
            # same coords. Escalate as drifts accumulate:
            #   tier 1 (drift 1-2): 50 ps @ 1.0 fs, taup=5   (slow-drift fix)
            #   tier 2 (drift 3-4): 100 ps @ 0.5 fs, taup=2, restart from TRUE eq
            #                       (finer integration + strong thermostat walks
            #                        clashing atoms apart from a pristine start)
            #   tier 3 (drift >=5): CPU density-settle (robust, integrator-agnostic)
            # ----------------------------------------------------------------
            # v2.5.60 FIX: distinguish a GENUINE STERIC CLASH from a BENIGN box
            # fluctuation before choosing the restart source.
            #
            # pmemd.cuda raises "box dimensions changed too much" whenever the
            # NPT box edge moves past its grid-cell tolerance -- AND THE AMBER
            # DOCS EXPLICITLY NOTE THIS FIRES EVEN WHEN THE SYSTEM IS HEALTHY.
            # The prescribed remedy is to restart FROM THE LATEST restart file
            # so a fresh grid is built at the CURRENT state (forward progress).
            #
            # The old ladder forced tier-2 to rewind to orig_eq.rst once
            # _drift_count>=3 and then FROZE there (drifts 3..10 all restarted
            # from the pristine origin), regenerating the same origin grid that
            # tripped the check -> 0 forward progress -> ladder exhausted. That
            # rewind is the RIGHT medicine for a real clash (lambda=0.150: VDW
            # spike + PRESS>8000) but the WRONG medicine for a benign
            # fluctuation (density ~1.01, no spike), which just needs to advance.
            #
            # Classifier (reads the prod.out we just halted on):
            #   clash  = last Density outside [0.85,1.15]  OR  PRESS magnitude
            #            > 5000  OR  a VDWAALS field that is NaN/overflow.
            #   benign = otherwise (box fine, just past grid tolerance).
            _last_dens_d=$(grep -oE 'Density[[:space:]]*=[[:space:]]*[0-9.]+' prod.out 2>/dev/null | tail -1 | grep -oE '[0-9.]+$') || true
            _last_press=$(grep -oE 'PRESS =[[:space:]]*[-0-9.]+' prod.out 2>/dev/null | tail -1 | grep -oE '[-0-9.]+$') || true
            _is_clash=0
            if grep -qiE 'VDWAALS[[:space:]]*=[[:space:]]*(\*+|[^ ]*[Nn]a[Nn])' prod.out 2>/dev/null; then _is_clash=1; fi
            if [ -n "${{_last_dens_d:-}}" ]; then
                _dok=$(awk -v d="$_last_dens_d" 'BEGIN{{print (d>=0.85 && d<=1.15)?1:0}}')
                if [ "$_dok" -ne 1 ]; then _is_clash=1; fi
            fi
            if [ -n "${{_last_press:-}}" ]; then
                _pbad=$(awk -v p="$_last_press" 'BEGIN{{ap=(p<0)?-p:p; print (ap>5000)?1:0}}')
                if [ "$_pbad" -eq 1 ]; then _is_clash=1; fi
            fi

            # Defaults = tier-1 gentle regrid from the LATEST good restart.
            local _rg_nst=50000 _rg_dt=0.001 _rg_taup=5.0 _rg_src="$good_rst" _tier=1
            if [ "$_is_clash" -eq 1 ]; then
                # Genuine clash: escalate to the fine-integration tier and, on
                # the FIRST clash escalation only, rewind to the pristine origin
                # to walk the clashing atoms apart. Subsequent clash attempts
                # resume from the latest good restart (the settled checkpoint),
                # so we never freeze re-rewinding to the same origin.
                if [ "$_drift_count" -ge 3 ]; then
                    _rg_nst=100000; _rg_dt=0.0005; _rg_taup=2.0; _tier=2
                    if [ "${{_clash_rewound:-0}}" -eq 0 ]; then
                        _rg_src="orig_eq.rst"; _clash_rewound=1
                    else
                        _rg_src="$good_rst"
                    fi
                fi
            else
                # Benign fluctuation: NEVER rewind to origin. Keep restarting
                # from the latest good restart so each regrid advances the
                # trajectory (Amber's documented remedy). Tighten integration a
                # little once drifts persist, but stay forward-progressing.
                if [ "$_drift_count" -ge 3 ]; then
                    _rg_nst=100000; _rg_dt=0.0005; _rg_taup=2.0; _tier=2
                fi
                _rg_src="$good_rst"
            fi
            _rec "[run_stage] lambda=${{LAM}} drift class=$( [ "$_is_clash" -eq 1 ] && echo CLASH || echo BENIGN ) (density=${{_last_dens_d:-?}}, press=${{_last_press:-?}})"
            _rec "[run_stage] lambda=${{LAM}} box drift #${{_drift_count}} on attempt ${{attempt}}; GPU regrid tier=${{_tier}} (${{_rg_nst}} steps @ ${{_rg_dt}} ps) from ${{_rg_src}}"
            if [ "$_drift_count" -ge 5 ] && [ "$cpu_used" -eq 0 ]; then
                # tier 3: escalate to CPU density-settle before the next GPU try
                cpu_density_settle || true
            fi
            # clean-template regrid input (no fragile sed on prod.in)
            _write_recovery_mdin prod_regrid.in "$_rg_nst" "$_rg_dt" "$_rg_taup" "GPU regrid tier ${{_tier}}"
            set +e
            pmemd.cuda -O -i prod_regrid.in -p system.prmtop -c "$_rg_src" \
                       -o prod_regrid.out.${{attempt}} -r prod_regrid.rst -x prod_regrid.nc \
                       > prod_regrid.console.${{attempt}} 2>&1
            regrid_rc=$?; set -e
            if [ "$regrid_rc" -eq 0 ] && [ -s prod_regrid.rst ]; then
                cp -f prod_regrid.rst "$good_rst"; rm -f prod.rst
                _rec "[run_stage] lambda=${{LAM}} GPU regrid OK (tier ${{_tier}}); resuming prod"
            else
                _rec "[run_stage] lambda=${{LAM}} GPU regrid FAILED (tier ${{_tier}}, rc=${{regrid_rc}})"
                tail -n 30 prod_regrid.console.${{attempt}} >&2 || true
                if [ "$attempt" -ge 2 ] && [ "$cpu_used" -eq 0 ]; then
                    cpu_density_settle || true
                fi
            fi
            attempt=$((attempt+1)); continue
        fi
        echo "============================================================" >&2
        _rec "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=prod rc=${{rc}} (non-box error)"
        grep -iE 'STOP PMEMD|Terminated Abnormally|NaN|vlimit exceeded|Input errors|ERROR' \
             prod.console.${{attempt}} prod.out 2>/dev/null | tail -n 20 >&2 || true
        tail -n 15 prod.out >&2 || true
        echo "============================================================" >&2
        return "$rc"
    done
    if [ "$cpu_used" -eq 1 ]; then
        _rec "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=prod -- still drifting after ${{max}} restarts INCLUDING a CPU density-settle"
    else
        _rec "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=prod -- box still drifting after ${{max}} restarts"
    fi
    return 255
}}
run_prod_with_restart
echo "[run_stage] lambda=${{LAM}} all stages completed OK"
"""
        script.write_text(h + body); script.chmod(0o755)
        return script

    def build_lsf_hremd(self, leg_dir, leg_name):
        nwin = len(self._active_lambdas)   # v2.5.73: was cfg.lambdas (wrong count)
        script = leg_dir / f"run_{leg_name}.lsf"
        h  = "#!/bin/bash\n"
        h += lsf_banner(f"{leg_name} HREMD")
        h += f"#BSUB -q {self.hpc_cfg.queue_gpu}\n"
        h += f"#BSUB -P {self.hpc_cfg.project}\n"
        h += f"#BSUB -J fep_{leg_name}_hremd\n"
        h += f"#BSUB -W {self.hpc_cfg.walltime}\n"
        h += f"#BSUB -o fep_{leg_name}.hremd.%J.out\n"
        h += f"#BSUB -e fep_{leg_name}.hremd.%J.err\n"
        _mem = getattr(self.hpc_cfg, "fep_mem_mb", 8192)
        h += f"#BSUB -n {nwin}\n"
        h += '#BSUB -R "span[hosts=1]"\n'
        h += f'#BSUB -R "rusage[mem={_mem}]"\n'
        h += self._avoid_select_line()
        h += '#BSUB -gpu "num=1"\n'
        h += '#BSUB -R "span[ptile=' + str(min(nwin, 4)) + ']"\n\n'
        h += "module purge\n"
        for m in self.hpc_cfg.modules: h += f"module load {m}\n"
        if self.hpc_cfg.venv_activate: h += f"source {self.hpc_cfg.venv_activate}\n"
        h += "set -euo pipefail\n\n"
        groupfile_path = leg_dir / "groupfile.prod"
        gf_lines = []
        for lam in self._active_lambdas:   # v2.5.73: was cfg.lambdas
            wd = f"{leg_dir}/lambda_{lam:.3f}"
            gf_lines.append(
                f"-O -i {wd}/prod.in -p {wd}/system.prmtop -c {wd}/eq.rst "
                f"-o {wd}/prod.out -r {wd}/prod.rst -x {wd}/prod.nc "
                f"-inf {wd}/prod.mdinfo")
        groupfile_path.write_text("\n".join(gf_lines) + "\n")
        lam_array = " ".join(f"{l:.3f}" for l in self._active_lambdas)   # v2.5.73
        body = f"""LAMBDAS=({lam_array})

# --- v2.4.15: surface the REAL pmemd failure instead of a bare exit 1 ---
run_stage() {{
    local stage="$1"; shift
    local mdout="" prev=""
    for a in "$@"; do
        if [ "$prev" = "-o" ]; then mdout="$a"; fi
        prev="$a"
    done
    echo "[run_stage] hremd ${{stage}} starting: $*"
    set +e; "$@"; local rc=$?; set -e
    if [ "$rc" -ne 0 ]; then
        echo "============================================================" >&2
        echo "STAGE FAILED: leg={leg_name} (hremd) stage=${{stage}} rc=${{rc}}" >&2
        if [ -n "$mdout" ] && [ -f "$mdout" ]; then
            echo "---- cause (from ${{mdout}}) ----" >&2
            grep -iE 'STOP PMEMD|Terminated Abnormally|NaN|vlimit exceeded|Input errors|ERROR' \\
                 "$mdout" | tail -n 20 >&2 || true
            echo "---- last 15 lines of ${{mdout}} ----" >&2
            tail -n 15 "$mdout" >&2 || true
        else
            echo "(no mdout '$mdout' produced -- pmemd died before output)" >&2
        fi
        echo "============================================================" >&2
        exit "$rc"
    fi
}}

# Per-window minimization + equilibration before replica exchange.
for LAM in "${{LAMBDAS[@]}}"; do
  WD={leg_dir}/lambda_${{LAM}}
  cd "$WD"
  run_stage "min(lam=${{LAM}})" pmemd.cuda -O -i min.in -p system.prmtop -c system.inpcrd \\
            -o min.out -r min.rst -x min.nc
  if [ -f heat.in ]; then
  run_stage "heat(lam=${{LAM}})" pmemd.cuda -O -i heat.in -p system.prmtop -c min.rst -ref min.rst \\
            -o heat.out -r heat.rst -x heat.nc
  HEAT_TEMP_MAX_K=${{HEAT_TEMP_MAX_K:-1000}}
  if [ -f heat.out ]; then
      maxTh=$(set +o pipefail; set +e; \\
             grep -aoE 'TEMP\\(K\\) *= *[-0-9.]+' heat.out 2>/dev/null \\
             | grep -aoE '[-0-9.]+$' \\
             | awk 'BEGIN{{m=""}} {{if(m==""||$1>m)m=$1}} END{{if(m!="")print m}}'; \\
             true)
      if [ -n "$maxTh" ] && awk -v t="$maxTh" -v m="$HEAT_TEMP_MAX_K" 'BEGIN{{exit !(t>m)}}'; then
          echo "STAGE FAILED: leg={leg_name} lambda=${{LAM}} stage=heat INSTABILITY: peak TEMP=${{maxTh}} K > ${{HEAT_TEMP_MAX_K}} K" >&2
          exit 70
      fi
  fi
  DENS_C=heat.rst
  else
  DENS_C=min.rst
  fi
  run_stage "dens(lam=${{LAM}})" pmemd.cuda -O -i dens.in -p system.prmtop -c "$DENS_C" -ref "$DENS_C" \\
            -o dens.out -r dens.rst -x dens.nc
  run_stage "eq(lam=${{LAM}})"  pmemd.cuda -O -i eq.in  -p system.prmtop -c dens.rst -ref dens.rst \\
            -o eq.out  -r eq.rst  -x eq.nc
  # --- final60 BUG 1 FIX: per-window equilibration stability gate (HREMD) ---
  # Same gate as the LSF-array path: a soft-core blow-up (e.g. complex_vdw
  # lambda~0.50 -> T=14,963 K) exits 0 but writes a poisoned eq.rst that would
  # then corrupt every replica in the exchange. Fail this window if eq peaked
  # above EQ_TEMP_MAX_K so it never enters replica exchange.
  EQ_TEMP_MAX_K=${{EQ_TEMP_MAX_K:-1000}}
  if [ -f eq.out ]; then
      # v2.5.7 FIX: SIGPIPE-safe peak-temp scan (see array builder note above).
      maxT=$(set +o pipefail; set +e; \
             grep -aoE 'TEMP\\(K\\) *= *[-0-9.]+' eq.out 2>/dev/null \
             | grep -aoE '[-0-9.]+$' \
             | awk 'BEGIN{{m=""}} {{if(m==""||$1>m)m=$1}} END{{if(m!="")print m}}'; \
             true)
      if [ -n "$maxT" ] && awk -v t="$maxT" -v m="$EQ_TEMP_MAX_K" 'BEGIN{{exit !(t>m)}}'; then
          echo "============================================================" >&2
          echo "STAGE FAILED: leg={leg_name} (hremd) lambda=${{LAM}} stage=eq INSTABILITY: peak TEMP=${{maxT}} K > ${{EQ_TEMP_MAX_K}} K" >&2
          grep -aiE 'vlimit exceeded|NaN|Terminated Abnormally' eq.out | tail -n 20 >&2 || true
          echo "(eq integrated to a non-physical temperature; refusing to enter" >&2
          echo " replica exchange on a blown-up window.)" >&2
          echo "============================================================" >&2
          exit 71
      fi
  fi
done

# Replica-exchange production (guarded).
cd {leg_dir}
echo "[run_stage] hremd mpirun replica-exchange starting ({nwin} replicas)"
set +e
mpirun -np {nwin} pmemd.cuda.MPI -ng {nwin} \\
       -groupfile {groupfile_path} -rem 3 -remlog rem.log
rc=$?
set -e
if [ "$rc" -ne 0 ]; then
    echo "============================================================" >&2
    echo "HREMD PRODUCTION FAILED: leg={leg_name} rc=${{rc}}" >&2
    if [ -f rem.log ]; then
        echo "---- last 20 lines of rem.log ----" >&2
        tail -n 20 rem.log >&2 || true
    fi
    echo "---- per-window prod.out crash scan ----" >&2
    for LAM in "${{LAMBDAS[@]}}"; do
        po="{leg_dir}/lambda_${{LAM}}/prod.out"
        if [ -f "$po" ]; then
            hit=$(grep -iE 'STOP PMEMD|Terminated Abnormally|NaN|vlimit exceeded|ERROR' "$po" | tail -n 3 || true)
            if [ -n "$hit" ]; then echo "  lambda=${{LAM}}: $hit" >&2; fi
        else
            echo "  lambda=${{LAM}}: no prod.out produced" >&2
        fi
    done
    echo "============================================================" >&2
    exit "$rc"
fi
echo "[run_stage] hremd all stages completed OK"
"""
        script.write_text(h + body); script.chmod(0o755)
        return script

    # PATCH A: analyzer per leg
    def _build_analyze_lsf(self, leg_dir, leg_name):
        from .submit import _WORKFLOW_ROOT
        queue = getattr(self.hpc_cfg, "queue_cpu", self.hpc_cfg.queue_gpu)
        lam_csv = ",".join(f"{l:.6f}" for l in self._active_lambdas)
        script = leg_dir / f"analyze_{leg_name}.lsf"
        h  = "#!/bin/bash\n"
        h += lsf_banner(f"analyze {leg_name}")
        h += f"#BSUB -q {queue}\n"
        h += f"#BSUB -P {self.hpc_cfg.project}\n"
        h += f"#BSUB -J fep_{leg_name}_analyze\n"
        h += "#BSUB -W 04:00\n#BSUB -n 1\n"
        h += f"#BSUB -o analyze_{leg_name}.%J.out\n"
        h += f"#BSUB -e analyze_{leg_name}.%J.err\n\n"
        h += "module purge\n"
        for m in self.hpc_cfg.modules: h += f"module load {m}\n"
        if self.hpc_cfg.venv_activate: h += f"source {self.hpc_cfg.venv_activate}\n"
        h += f'export PYTHONPATH="{_WORKFLOW_ROOT}${{PYTHONPATH:+:$PYTHONPATH}}"\n'
        h += "set -uo pipefail\n\n"
        temp = self.cfg.temperature_K
        body = f"""python - <<'PY'
import json, sys
from pathlib import Path
from amber_md.fep import FEPAnalyzer, analyze_restraint_leg_optionA
leg_dir = Path("{leg_dir}")
lambdas = [float(x) for x in "{lam_csv}".split(",")]
try:
    # v2.5.63: Option-A restraint leg (plain MD, no clambda) scored
    # analytically from the Boresch correction; MBAR only for TI legs.
    res = analyze_restraint_leg_optionA(
        leg_dir, lambdas, temperature_K={self.cfg.temperature_K})
    if res is None:
        res = FEPAnalyzer(leg_dir, lambdas,
                          temperature_K={self.cfg.temperature_K}).run()
except Exception as e:
    import traceback; traceback.print_exc()
    print("ANALYZER FAILED: {leg_name}: %s" % e)
    sys.exit(3)

(leg_dir / "summary.json").write_text(json.dumps(res, indent=2, default=str))
print("======== {leg_name} leg ========")
print("headline estimator :", res.get("headline_estimator"))
print("estimator used     :", res.get("estimator_used"))
print("dG (kcal/mol)      :", res.get("dG_kcal_mol"))
print("windows            : %s/%s" % (res.get("n_windows"), res.get("n_requested")))
if res.get("missing_windows"):
    print("missing windows    :", res.get("missing_windows"))
if res.get("dG_boresch_correction") is not None:
    print("  uncorrected dG  :", res.get("dG_uncorrected_kcal_mol"))
    print("  Boresch correction:", res.get("dG_boresch_correction"))

if res.get("dG_kcal_mol") is None:
    print("ANALYZER STATUS: FAILED (no dG produced)")
    sys.exit(2)
if not res.get("complete", True):
    print("ANALYZER STATUS: INCOMPLETE (%s/%s windows; result UNRELIABLE)"
          % (res.get("n_windows"), res.get("n_requested")))
    sys.exit(1)
print("ANALYZER STATUS: OK")
PY
rc=$?
echo "[analyze] {leg_name} python exit code: $rc"

# v2.5.73: also generate the TI/BAR/MBAR convergence CSV+PNG (non-fatal).
_conv="{_WORKFLOW_ROOT}/tools/convergence_analysis.py"
if [ -f "$_conv" ]; then
  echo "[analyze] {leg_name} running convergence_analysis.py ..."
  python "$_conv" "{leg_dir}" --temp {temp} --estimator auto \
      --out "{leg_dir}/convergence" || echo "[analyze] {leg_name} convergence failed (non-fatal)"
else
  echo "[analyze] {leg_name} WARN: $_conv not found; skipped convergence CSV"
fi

exit $rc
"""
        script.write_text(h + body); script.chmod(0o755)
        return script

    # v2.4.14 NEW: cycle-closer
    def build_cycle_closer_lsf(self, fep_root):
        """Single-CPU LSF job that waits for BOTH analyzers, then writes
        fep/ABFE_RESULT.txt + ABFE_RESULT.json with cycle-closed ΔG_bind.
        """
        from .submit import _WORKFLOW_ROOT
        fep_root = Path(fep_root)
        queue = getattr(self.hpc_cfg, "queue_cpu", self.hpc_cfg.queue_gpu)
        script = fep_root / "cycle_close_abfe.lsf"
        h  = "#!/bin/bash\n"
        h += lsf_banner("ABFE cycle-closer")
        h += f"#BSUB -q {queue}\n"
        h += f"#BSUB -P {self.hpc_cfg.project}\n"
        h += "#BSUB -J fep_cycle_close\n"
        h += "#BSUB -W 00:30\n#BSUB -n 1\n"
        h += "#BSUB -o cycle_close.%J.out\n"
        h += "#BSUB -e cycle_close.%J.err\n\n"
        h += "module purge\n"
        for m in self.hpc_cfg.modules: h += f"module load {m}\n"
        if self.hpc_cfg.venv_activate: h += f"source {self.hpc_cfg.venv_activate}\n"
        h += f'export PYTHONPATH="{_WORKFLOW_ROOT}${{PYTHONPATH:+:$PYTHONPATH}}"\n'
        h += "set -uo pipefail\n\n"
        body = f"""python - <<'PY'
import json
from pathlib import Path
from amber_md.fep import FEPAnalyzer

import sys
root = Path("{fep_root}")
T    = {self.cfg.temperature_K}
DECHG = [{', '.join(f'{l:.6f}' for l in getattr(self.cfg, "decharge_lambdas", self.cfg.lambdas))}]
VDW   = [{', '.join(f'{l:.6f}' for l in getattr(self.cfg, "vdw_lambdas",      self.cfg.lambdas))}]
REST  = [{', '.join(f'{l:.6f}' for l in (self.cfg.restraint_lambdas_fine if getattr(self.cfg,"use_fine_restraint_lambdas",False) else getattr(self.cfg,"restraint_lambdas",self.cfg.lambdas)))}]

def _run(sub, lams):
    d = root / sub
    if not d.exists(): return None
    s = d / "summary.json"
    # PATCH Bug 4 (v2.4.19): only trust a cached summary if it is COMPLETE;
    # otherwise re-analyze (resume may have added windows since it was written).
    if s.exists():
        try:
            cached = json.loads(s.read_text())
            if cached.get("complete", False):
                return cached
        except Exception:
            pass
    res = FEPAnalyzer(d, lams, temperature_K=T).run()
    s.write_text(json.dumps(res, indent=2, default=str))
    return res

def _leg(prefix, with_restraint=False):
    # v2.5.16: complex leg = decharge + vdw + restraint(MD) + Boresch(analytic).
    # v2.5.38: also RETURN the restraint-leg result so the cycle-closer can
    # enforce its reliability gate (dG_reliable) on the headline dG_bind.
    dchg = _run(prefix + "_decharge", DECHG)
    vdw  = _run(prefix + "_vdw",      VDW)
    if dchg is None or vdw is None: return None, dchg, vdw, None
    a = dchg.get("dG_kcal_mol"); b = vdw.get("dG_kcal_mol")
    if a is None or b is None: return None, dchg, vdw, None
    total = a + b
    rest = None
    if with_restraint:
        rest = _run(prefix + "_restraint", REST)
        if rest is None or rest.get("dG_kcal_mol") is None:
            return None, dchg, vdw, rest
        total += rest.get("dG_kcal_mol")          # MD restraint free energy
        bcorr = rest.get("dG_boresch_correction")
        if bcorr is not None: total += bcorr
    return total, dchg, vdw, rest

# v2.5.18: AUTOMATIC Rocklin PB charge correction for net-charged ligands.
# Computes the complex-minus-solvent finite-size correction from the final
# frame of each decoupled-end-state leg and writes charge_correction.json,
# which is then consumed below. Self-skips (0.0) for neutral ligands or if
# any dependency / input is unavailable -- never aborts the cycle-closer.
_ccf = root / "charge_correction.json"
if not _ccf.exists():
    try:
        import parmed as _pmd
        _cprm = {getattr(self.cfg, "complex_prmtop", "") or ""!r}
        _sprm = {getattr(self.cfg, "solvent_prmtop", "") or ""!r}
        _lrn  = {getattr(self.cfg, "ligand_resname", "LIG")!r}
        _smask = {getattr(self.cfg, "solvent_mask", ":WAT,K+,Na+,Cl-")!r}
        _water = {getattr(self.cfg, "water_model", "tip3p")!r}
        _T = T
        _eps = {{"tip3p": 97.0, "spce": 97.0, "tip4pew": 97.0}}.get(str(_water).lower(), 97.0)
        _qlig = 0.0
        if _cprm and Path(_cprm).exists():
            _m = _pmd.load_file(str(_cprm))
            _lr = [r for r in _m.residues if r.name.strip() == _lrn]
            if _lr:
                _qlig = round(sum(a.charge for a in _lr[0].atoms))
        if abs(_qlig) >= 1:
            from amber_md.charge_correction import (compute_charge_correction,
                                                    align_complex)
            def _final_crd(sub):
                d = root / sub
                for cand in (sorted(d.glob("lambda_*/prod.rst7")) +
                             sorted(d.glob("lambda_*/prod.rst")) +
                             sorted(d.glob("lambda_*/prod.nc"))):
                    if cand.exists():
                        return cand
                return None
            _ccomp = _final_crd("complex_decharge")
            _csolv = _final_crd("solvent_decharge")
            _dG_c = _dG_s = 0.0
            if _cprm and _ccomp:
                _a = align_complex(str(_cprm), str(_ccomp), _smask)
                _dfc = compute_charge_correction(str(_cprm), _a, _lrn, _smask, _T, _eps)
                _dG_c = float(_dfc["Total"].iloc[-1]) if hasattr(_dfc, "columns") else float(_dfc)
            if _sprm and _csolv:
                _a = align_complex(str(_sprm), str(_csolv), _smask)
                _dfs = compute_charge_correction(str(_sprm), _a, _lrn, _smask, _T, _eps)
                _dG_s = float(_dfs["Total"].iloc[-1]) if hasattr(_dfs, "columns") else float(_dfs)
            _ccf.write_text(json.dumps(
                {{"dG_charge_correction_kcal_mol": _dG_c - _dG_s,
                  "complex_total": _dG_c, "solvent_total": _dG_s,
                  "ligand_net_charge": _qlig}}, indent=2))
            print("[charge-correction] q=%+d -> dG_corr=%+.3f kcal/mol"
                  % (int(_qlig), _dG_c - _dG_s))
        else:
            print("[charge-correction] neutral ligand; no correction needed.")
    except Exception as _e:
        print("[charge-correction] skipped (%s); using 0.0" % _e)
dG_charge = 0.0
if _ccf.exists():
    try: dG_charge = float(json.loads(_ccf.read_text()).get("dG_charge_correction_kcal_mol", 0.0))
    except Exception: dG_charge = 0.0

ctot, cdc, cvd, crest = _leg("complex", with_restraint=True)
stot, sdc, svd, _srest = _leg("solvent")
dG_bind = (None if (ctot is None or stot is None)
           else -(ctot - stot) + dG_charge)

def Fmt(x): return "(failed)" if x is None else ("%+8.3f kcal/mol" % x)
g = lambda r: None if r is None else r.get("dG_kcal_mol")

# v2.4.15: assess completeness of every contributing leg.
_subres = {{"complex_decharge": cdc, "complex_vdw": cvd,
            "complex_restraint": crest,
            "solvent_decharge": sdc, "solvent_vdw": svd}}
def _leg_ok(r):
    return (r is not None and r.get("dG_kcal_mol") is not None
            and r.get("complete", True))
incomplete = sorted(k for k, r in _subres.items() if not _leg_ok(r))
_DG_BIND_SANITY_KCAL = 25.0
nonconv = sorted(k for k, r in _subres.items()
                 if r is not None and r.get("estimator_converged") is False)
# final40 FIX 2: a leg where BAR *and* MBAR failed/diverged (only TI survived)
# is NOT trustworthy even though its (TI-only) spread is small. Track these
# separately so the cycle is honestly untrusted with the REAL reason instead
# of a phantom huge spread or a false GO.
mbar_failed = sorted(k for k, r in _subres.items()
                     if r is not None and r.get("only_ti_survived"))
sane_mag = (dG_bind is not None) and (abs(dG_bind) <= _DG_BIND_SANITY_KCAL)
# v2.5.38: a leg flagged dG_reliable=False by the analyzer reliability gate
# (max|dV/dl|>200 or |BAR-TI|>50) must untrust the whole cycle.
unreliable = sorted(k for k, r in _subres.items()
                    if r is not None and r.get("dG_reliable") is False)
all_ok = ((dG_bind is not None) and (len(incomplete) == 0)
          and (len(nonconv) == 0) and (len(mbar_failed) == 0)
          and (len(unreliable) == 0) and sane_mag)

def Cmpl(r):
    if r is None: return "MISSING"
    if r.get("dG_kcal_mol") is None: return "FAILED"
    if not r.get("complete", True):
        return "INCOMPLETE %s/%s" % (r.get("n_windows"), r.get("n_requested"))
    return "complete"

lines = ["============== ABFE RESULT (two-stage) =============="]
lines += ["  COMPLEX:  decharge=" + Fmt(g(cdc)) + "  vdw=" + Fmt(g(cvd)) + "  total(+Boresch)=" + Fmt(ctot)]
lines += ["  SOLVENT:  decharge=" + Fmt(g(sdc)) + "  vdw=" + Fmt(g(svd)) + "  total=" + Fmt(stot)]
lines += ["  leg status: complex_decharge=" + Cmpl(cdc) + ", complex_vdw=" + Cmpl(cvd)]
lines += ["              complex_restraint=" + Cmpl(crest) +
          ("" if not crest or crest.get("dG_reliable", True)
           else " [UNRELIABLE: " + "; ".join(crest.get("reliability_reasons", []) or ["gate failed"]) + "]")]
lines += ["              solvent_decharge=" + Cmpl(sdc) + ", solvent_vdw=" + Cmpl(svd)]
lines += ["-----------------------------------------------------"]
if all_ok:
    lines += ["  dG_bind = " + Fmt(dG_bind)]
else:
    why_bits = []
    if dG_bind is None:        why_bits.append("missing/failed legs")
    if incomplete:             why_bits.append("incomplete: " + ", ".join(incomplete))
    if nonconv:                why_bits.append("TI/BAR/MBAR disagree: " + ", ".join(nonconv))
    if dG_bind is not None and abs(dG_bind) > _DG_BIND_SANITY_KCAL:
        why_bits.append("dG_bind magnitude %.1f > %.0f (unphysical)"
                        % (abs(dG_bind), _DG_BIND_SANITY_KCAL))
    if mbar_failed:
        why_bits.append("MBAR/BAR failed on %s (TI-only; cross-energy "
                        "estimator unavailable -- typical of truncated / "
                        "too-few-window runs)" % ", ".join(mbar_failed))
    if unreliable:
        why_bits.append("RELIABILITY GATE failed on %s (runaway dV/dl or "
                        "BAR/TI disagreement -- physically broken leg)"
                        % ", ".join(unreliable))
    why = "; ".join(why_bits) if why_bits else "untrusted"
    lines += ["  dG_bind = " + Fmt(dG_bind) + "   *** UNTRUSTED (" + why + ") ***"]
lines += ["====================================================="]
out = "\\n".join(lines) + "\\n"
print(out)
(root / "ABFE_RESULT.txt").write_text(out)
(root / "ABFE_RESULT.json").write_text(json.dumps({{
    "dG_complex_plus_restr_kcal_mol": ctot,
    "dG_solvent_kcal_mol": stot,
    "dG_bind_kcal_mol": dG_bind,
    "trusted": all_ok,
    "incomplete_legs": incomplete,
    "mbar_failed_legs": mbar_failed,
    "unreliable_legs": unreliable,
    "complex_restraint": g(crest),
    "restraint_reliable": (None if crest is None else crest.get("dG_reliable")),
    "complex_decharge": g(cdc), "complex_vdw": g(cvd),
    "solvent_decharge": g(sdc), "solvent_vdw": g(svd),
    "absolute_headline": None if cvd is None else cvd.get("headline_estimator"),
    "solvent_headline":  None if svd is None else svd.get("headline_estimator"),
    "estimator_spread_kcal": {{k: (None if r is None else r.get("estimator_spread_kcal"))
                               for k, r in _subres.items()}},
    "nonconverged_legs": nonconv,
    "temperature_K": T,
}}, indent=2))
if not all_ok:
    print("CYCLE-CLOSER STATUS: UNTRUSTED -- see leg status above")
    sys.exit(1)
print("CYCLE-CLOSER STATUS: OK")
PY
"""
        script.write_text(h + body); script.chmod(0o755)
        return script

    # ---------------------------------------------------------------
    # v2.4.18 NEW: resume support -- detect & resubmit only the windows
    # whose production did NOT finish (e.g. killed by HPC maintenance or
    # preemption). A window is "complete" iff its prod.out exists and
    # contains pmemd's end-of-run marker.
    # v2.4.18a: corrected to match the real _bsub_submit signature
    #   _bsub_submit(script, cwd, *, project, queue, walltime,
    #                extra_args=None) -> bare job-id str
    # (walltime is REQUIRED; dependencies go through extra_args=["-w",...]).
    # ---------------------------------------------------------------
    _PROD_DONE_MARKERS = ("Final Performance Info", "TIMINGS",
                          "Total wall time")

    def _window_complete(self, leg_dir, lam):
        po = Path(leg_dir) / f"lambda_{lam:.3f}" / "prod.out"
        if not po.exists() or po.stat().st_size == 0:
            return False
        try:
            txt = po.read_text(errors="replace")
        except Exception:
            return False
        return any(m in txt for m in self._PROD_DONE_MARKERS)

    def incomplete_indices(self, leg_dir):
        """1-based LSF array indices of windows that are NOT complete.
        Index i maps to self._active_lambdas[i-1], matching the
        LAMBDAS=(...) / LSB_JOBINDEX scheme in build_lsf_array()."""
        idx = []
        for i, lam in enumerate(self._active_lambdas, start=1):
            if not self._window_complete(leg_dir, lam):
                idx.append(i)
        return idx

    @staticmethod
    def _indices_to_bsub_ranges(indices):
        """Compress [2,3,4,7,9,10] -> '2-4,7,9-10' for #BSUB -J name[...]."""
        if not indices:
            return ""
        indices = sorted(set(indices))
        parts, start, prev = [], indices[0], indices[0]
        for n in indices[1:]:
            if n == prev + 1:
                prev = n; continue
            parts.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = prev = n
        parts.append(f"{start}-{prev}" if start != prev else f"{start}")
        return ",".join(parts)

    def build_lsf_array_resume(self, leg_dir, leg_name, indices):
        """Like build_lsf_array but the -J array spans only `indices`
        (the unfinished windows). Body is identical, so finished windows
        are simply never launched."""
        full = self.build_lsf_array(leg_dir, leg_name)
        txt = full.read_text()
        rng = self._indices_to_bsub_ranges(indices)
        txt = re.sub(r"#BSUB -J fep_" + re.escape(leg_name) + r"\[[^\]]*\]",
                     f"#BSUB -J fep_{leg_name}[{rng}]", txt, count=1)
        script = leg_dir / f"run_{leg_name}_resume.lsf"
        script.write_text(txt)
        script.chmod(0o755)
        return script

    def submit_leg_resume(self, leg_dir, leg_name):
        """Resubmit only unfinished windows + (re)submit the analyzer so
        the cycle-closer can fire once they complete. Returns the same
        dict shape as submit_leg: {'md':jid, 'analyze':jid}.

        Mirrors submit_leg()'s exact _bsub_submit usage (walltime
        required; dependency via extra_args=['-w','ended(JID)'])."""
        indices = self.incomplete_indices(leg_dir)
        md_jid = None
        if indices:
            rng = self._indices_to_bsub_ranges(indices)
            log.info("  resume %s: %d/%d windows incomplete -> array [%s]",
                     leg_name, len(indices), len(self._active_lambdas), rng)
            if self.hremd:
                # HREMD couples all replicas in one job; partial resume is
                # not meaningful -> full rerun.
                script = self.build_lsf_hremd(leg_dir, leg_name)
            else:
                script = self.build_lsf_array_resume(leg_dir, leg_name, indices)
            md_jid = _bsub_submit(
                script, leg_dir,
                project=self.hpc_cfg.project,
                queue=self.hpc_cfg.queue_gpu,
                walltime=self.hpc_cfg.walltime)
        else:
            log.info("  resume %s: all %d windows already complete; "
                     "re-running analyzer only.",
                     leg_name, len(self._active_lambdas))

        # Always (re)wire the analyzer so the cycle-closer can complete.
        # If we resubmitted MD, the analyzer waits on it; otherwise it can
        # run immediately against the already-finished windows.
        ana_script = self._build_analyze_lsf(leg_dir, leg_name)
        ana_queue = getattr(self.hpc_cfg, "queue_cpu", self.hpc_cfg.queue_gpu)
        extra = ["-w", f"ended({md_jid})"] if md_jid else None
        try:
            ana_jid = _bsub_submit(
                ana_script, leg_dir,
                project=self.hpc_cfg.project,
                queue=ana_queue, walltime="04:00",
                extra_args=extra)
        except Exception as e:
            log.warning("Could not submit analyzer for %s: %s", leg_name, e)
            ana_jid = None
        return {"md": md_jid, "analyze": ana_jid}

    def submit_leg(self, leg_dir, leg_name):
        script = (self.build_lsf_hremd(leg_dir, leg_name) if self.hremd
                  else self.build_lsf_array(leg_dir, leg_name))
        md_jid = _bsub_submit(
            script, leg_dir,
            project=self.hpc_cfg.project,
            queue=self.hpc_cfg.queue_gpu,
            walltime=self.hpc_cfg.walltime)
        ana_script = self._build_analyze_lsf(leg_dir, leg_name)
        ana_queue = getattr(self.hpc_cfg, "queue_cpu", self.hpc_cfg.queue_gpu)
        try:
            ana_jid = _bsub_submit(
                ana_script, leg_dir,
                project=self.hpc_cfg.project,
                queue=ana_queue, walltime="04:00",
                extra_args=["-w", f"ended({md_jid})"])
        except Exception as e:
            log.warning("Could not submit analyzer for %s: %s", leg_name, e)
            ana_jid = None
        return {"md": md_jid, "analyze": ana_jid}

    # v2.4.14 NEW
    def submit_cycle_closer(self, fep_root, analyzer_jids):
        # >>> PATCH Bug 1 (v2.4.19): read the TWO-STAGE analyzer JID keys
        #     (complex_decharge/complex_vdw/solvent_decharge/solvent_vdw),
        #     not the obsolete one-stage absolute/solvent keys. In two-stage
        #     mode the old keys were always absent -> guard tripped with
        #     "got abs=None sol=None" -> cycle-closer never submitted ->
        #     ABFE_RESULT.txt never written. Legacy one-stage keys kept as
        #     a fallback. Wait (ended()) on ALL FOUR per-stage analyzers.
        abs_jids = [analyzer_jids.get("complex_decharge"),
                    analyzer_jids.get("complex_vdw")]
        sol_jids = [analyzer_jids.get("solvent_decharge"),
                    analyzer_jids.get("solvent_vdw")]
        if not any(abs_jids):
            abs_jids = [analyzer_jids.get("absolute")]
        if not any(sol_jids):
            sol_jids = [analyzer_jids.get("solvent")]
        abs_jids = [j for j in abs_jids if j]
        sol_jids = [j for j in sol_jids if j]
        if not abs_jids or not sol_jids:
            log.warning("Cycle-closer not submitted: need both legs' "
                        "analyzer JIDs (got complex=%s solvent=%s)",
                        abs_jids, sol_jids)
            return None
        script = self.build_cycle_closer_lsf(fep_root)
        queue = getattr(self.hpc_cfg, "queue_cpu", self.hpc_cfg.queue_gpu)
        # ended() (not done()) so the closer still runs even when an analyzer
        # exits nonzero on an INCOMPLETE/FAILED leg; it then writes a clear
        # UNTRUSTED verdict instead of never firing.
        dep = " && ".join(f"ended({j})" for j in (abs_jids + sol_jids))
        try:
            jid = _bsub_submit(
                script, Path(fep_root),
                project=self.hpc_cfg.project,
                queue=queue, walltime="00:30",
                extra_args=["-w", dep])
            log.info("Cycle-closer queued JID=%s (waits [ended] on %s)",
                     jid, dep)
            return jid
        except Exception as e:
            log.warning("Could not submit cycle-closer: %s", e)
            return None
        # <<< PATCH Bug 1


# ----------------------------------------------------------------------


# === v2.5.63: Option-A restraint-leg analytic analyzer ==================
# The Option-A complex_restraint leg is plain NPT MD (prod.in has NO
# clambda/icfe), so alchemlyb/MBAR finds "no free energy section" and the
# leg yields dG=None -> analyze exited 2 (looked like a pipeline failure).
# Its free-energy contribution is the ANALYTIC Boresch correction, already
# written to boresch_correction.txt. This helper detects that mode and
# returns a COMPLETE summary using the analytic value, so analyze exits 0.
def analyze_restraint_leg_optionA(leg_dir, lambdas, temperature_K=298.0):
    """Summary dict for an Option-A (plain-MD) restraint leg.

    Returns None if the windows are a real TI restraint leg (icfe/clambda
    present) so the caller falls back to FEPAnalyzer/MBAR. Otherwise reports
    the analytic Boresch correction as the leg dG.
    """
    import re as _re
    from pathlib import Path
    leg_dir = Path(leg_dir)
    lams = [float(x) for x in lambdas]

    # Detect TI vs Option-A from produced prod.in/prod.out.
    is_ti, seen_any = False, False
    for lam in lams:
        d = leg_dir / ("lambda_%.3f" % lam)
        for fn in ("prod.in", "prod.out"):
            f = d / fn
            if f.exists():
                seen_any = True
                t = f.read_text(errors="ignore")
                if ("clambda" in t) or ("icfe=1" in t) or ("icfe =" in t):
                    is_ti = True
                break
    if is_ti or not seen_any:
        return None

    present, missing = [], []
    for lam in lams:
        po = leg_dir / ("lambda_%.3f" % lam) / "prod.out"
        ok = po.exists() and bool(
            _re.search(r"Final Performance Info|Total wall time",
                       po.read_text(errors="ignore")))
        (present if ok else missing).append(round(lam, 3))

    boresch = None
    cf = leg_dir / "boresch_correction.txt"
    if cf.exists():
        try: boresch = float(cf.read_text().strip())
        except ValueError: boresch = None

    complete = (len(missing) == 0) and (boresch is not None)
    return {
        "leg_type": "restraint_optionA",
        "method": "analytic_boresch",
        "dG_kcal_mol": (boresch if complete else None),
        "dG_uncorrected_kcal_mol": 0.0,
        "dG_boresch_correction": boresch,
        "estimators": {"analytic_boresch":
                       {"dG_kcal_mol": boresch, "err_kcal_mol": 0.0}},
        "headline_estimator": "analytic_boresch",
        "estimator_used": "analytic_boresch",
        "n_requested": len(lams),
        "n_windows": len(present),
        "missing_windows": missing,
        "complete": bool(complete),
    }
# === end v2.5.63 helper =================================================

class FEPAnalyzer:
    def __init__(self, leg_dir, lambdas, temperature_K=298.0,
                 gq_nodes=None, gq_weights=None):
        self.dir = Path(leg_dir); self.lambdas = list(lambdas); self.T = float(temperature_K)
        self._gq_nodes = list(gq_nodes) if gq_nodes else None
        self._gq_weights = list(gq_weights) if gq_weights else None

    # v2.5.48: FEPAnalyzer previously called self._window_complete(), but that
    # method (and its _PROD_DONE_MARKERS) lived ONLY on FEPSetup -- so the
    # analyzer crashed with AttributeError the first time a real prod.out
    # existed (lines _collect_dvdl / run). Define them here so FEPAnalyzer is
    # self-contained. A window counts as complete only if its prod.out carries
    # a genuine pmemd completion marker (not merely existing / mid-run).
    _PROD_DONE_MARKERS = ("Final Performance Info", "TIMINGS", "Total wall time")

    def _window_complete(self, leg_dir, lam):
        po = Path(leg_dir) / f"lambda_{lam:.3f}" / "prod.out"
        if not po.exists() or po.stat().st_size == 0:
            return False
        try:
            txt = po.read_text(errors="replace")
        except Exception:
            return False
        return any(m in txt for m in self._PROD_DONE_MARKERS)

    @staticmethod
    def _parse_dvdl(prod_out):
        """Return <DV/DL> from the AVERAGES block of a pmemd mdout.

        Hardened (v2.4.15): scans every DV/DL line within the averages
        section (not just the first), tolerates output with or without
        an "=" sign, and rejects non-finite values so a NaN/Inf window
        is reported as missing rather than silently integrated.
        """
        in_avg = False
        result = None
        for line in Path(prod_out).read_text().splitlines():
            if "A V E R A G E S   O V E R" in line:
                in_avg = True
                continue
            if not in_avg:
                continue
            if "R M S" in line:   # end of averages block
                break
            if "DV/DL" in line:
                tail = line.split("=", 1)[1] if "=" in line else line.split("DV/DL", 1)[1]
                for tok in tail.split():
                    try:
                        val = float(tok)
                    except ValueError:
                        continue
                    if math.isfinite(val):
                        result = val
                    break
        return result

    def _collect_dvdl(self):
        d = {}
        for lam in self.lambdas:
            out = self.dir / f"lambda_{lam:.3f}" / "prod.out"
            if not out.exists(): continue
            if not self._window_complete(self.dir, lam):
                continue
            v = self._parse_dvdl(out)
            if v is not None: d[lam] = v
        return d

    def _completeness(self, found):
        """Compare windows actually analyzed (found) against requested."""
        req = list(self.lambdas)
        fset = set(found)
        missing = [l for l in req if l not in fset]
        return {
            "n_requested":     len(req),
            "n_windows":       len(found),
            "missing_windows": [round(float(l), 3) for l in missing],
            "complete":        len(missing) == 0,
        }

    def _fallback_ti(self):
        d = self._collect_dvdl()
        comp = self._completeness(list(d.keys()))
        if not d:
            return {"dG_kcal_mol": None, "estimators": {}, "dvdl": {}, **comp}
        xs = sorted(d); ys = [d[x] for x in xs]
        try:
            import numpy as np
            dG = float(np.trapz(ys, xs))
        except ImportError:
            dG = sum(0.5*(ys[i]+ys[i+1])*(xs[i+1]-xs[i]) for i in range(len(xs)-1))
        if not comp["complete"]:
            log.warning("  INCOMPLETE leg %s: %d/%d windows; missing %s "
                        "(trapezoid TI over partial profile is UNRELIABLE)",
                        self.dir.name, comp["n_windows"], comp["n_requested"],
                        comp["missing_windows"])
        return {"dG_kcal_mol": dG,
                "estimators": {"TI_legacy": {"dG_kcal_mol": dG, "err_kcal_mol": None}},
                "dvdl": d, "headline_estimator": "TI_legacy", **comp}

    @staticmethod
    def _is_sane(e):
        if not e: return False
        dG = e.get("dG_kcal_mol")
        if dG is None: return False
        try:
            if not math.isfinite(dG) or abs(dG) > _HEADLINE_SANITY_MAX_KCAL: return False
        except (TypeError, ValueError): return False
        err = e.get("err_kcal_mol")
        if err is not None:
            try:
                if not math.isfinite(err): return False
            except (TypeError, ValueError): return False
        return True

    @staticmethod
    def _sanitize_u_nk(u_nk):
        """v2.5.14 FIX: CLIP-only u_nk sanitation -- drop NOTHING.

        Root cause (confirmed on real full-leg data, lig_12944901):
          * The ONLY corrupt cells are the decoupled end-state singularity:
            the fully-decoupled lambda evaluated from a coupled config prints
            either **** (Fortran field overflow -> non-finite) or a finite but
            physically-meaningless ~1e9 reduced energy. Across an entire leg
            this is a single column (the lambda=0.0 end state) in ~15-22% of
            frames; every other cell is clean (physical median ~O(1-100) kT).

        Why the v2.5.12/2.5.13 DROP strategies were wrong (both validated):
          * v2.5.12 dropped any FRAME with a bad cell -> emptied ~half the
            windows -> degenerate matrix -> MBAR dG=-710 vs TI=+8.9.
          * v2.5.13 dropped the bad COLUMN + residual rows -> on real data the
            1e4 ABSOLUTE cap flagged ~100% of cells (reduced energies are
            offset), nuking the whole matrix. Even calibrated, dropping a
            column removes a real thermodynamic state from the MBAR cycle.
          * Synthetic ground-truth test (12-state harmonic ladder, 22% cells
            replaced by inf / +1.6e9 / -8.3e6): clip-to-ceiling recovers the
            true dG to 0.03 kcal/mol while preserving ALL rows and columns;
            the old 1e4-cap+row-drop kept only 316/6000 frames and was 5x
            less accurate.

        Correct behaviour (Boltzmann limit of an infinite-energy state):
          A config with ~infinite reduced energy at some state contributes
          ZERO statistical weight there. exp(-u) with u >> 1 underflows to 0,
          which is exactly right. So we CLIP every non-finite cell, and every
          finite cell with |u_nk| > _UNK_SANE_MAX_KT, to +_UNK_SANE_MAX_KT
          (a large POSITIVE reduced energy => zero weight). The matrix stays
          square (k x k); no rows, columns, or sampled lambda groups are
          removed; TI is untouched.

        Returns (clean_u_nk, n_clipped_cells). Never raises (best-effort).
        """
        try:
            import numpy as np, pandas as pd
        except ImportError:
            return u_nk, 0
        vals = u_nk.to_numpy(dtype=float, copy=True)
        if vals.size == 0:
            return u_nk, 0
        finite = np.isfinite(vals)
        bad = ~finite
        bad |= np.where(finite, vals, np.inf) > _UNK_SANE_MAX_KT
        bad |= np.where(finite, vals, -np.inf) < -_UNK_SANE_MAX_KT
        n_clipped = int(bad.sum())
        if n_clipped == 0:
            return u_nk, 0
        vals[bad] = _UNK_SANE_MAX_KT
        clean = pd.DataFrame(vals, index=u_nk.index, columns=u_nk.columns)
        return clean, n_clipped
    @staticmethod
    def _u_nk_rank_cond(u_nk):
        """v2.5.15: numerical rank + condition number of u_nk (SVD on the
        column-centered matrix). Returns (rank, cond, n_states, well_posed).
        MBAR needs the column space to span all n_states; near-collinear
        soft-core columns make it under-determined (the -1477 artifact).
        BAR uses adjacent pairs only and is NOT gated by this. Best-effort:
        on any failure returns well_posed=True (legacy fallback).
        """
        try:
            import numpy as np
            v = u_nk.to_numpy(dtype=float, copy=True)
            n_states = v.shape[1]
            if v.size == 0 or n_states == 0:
                return 0, float('inf'), n_states, True
            vc = v - v.mean(axis=0, keepdims=True)
            sv = np.linalg.svd(vc, compute_uv=False)
            if sv[0] == 0:
                return 0, float('inf'), n_states, False
            rank = int((sv > _UNK_RANK_RTOL * sv[0]).sum())
            cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf')
            well = (rank >= n_states) and (cond <= _UNK_COND_MAX)
            return rank, cond, n_states, well
        except Exception:  # noqa: BLE001
            return -1, float('nan'), -1, True

    @staticmethod
    def _reconcile_u_nk_states(u_nk):
        """final60 BUG 2 FIX (guard): make the reduced-potential matrix square.

        Even with clamping, a window can legitimately contribute ZERO surviving
        frames (e.g. the alchemlyb parser produced no usable rows for one
        lambda). MBAR/BAR require the set of sampled lambda groups (index) to
        match the set of evaluated states (columns). When they differ we
        restrict BOTH to their intersection so the fit sees a (k x k) matrix
        instead of a rectangular one. Returns (u_nk, dropped_states) where
        dropped_states lists the column labels removed (for logging).
        """
        try:
            import numpy as np  # noqa: F401
        except ImportError:
            return u_nk, []
        try:
            sampled = list(dict.fromkeys(
                u_nk.index.get_level_values("lambdas")))
        except Exception:
            return u_nk, []
        cols = list(u_nk.columns)
        # column labels may be scalars or 1-tuples; normalise for comparison.
        def _key(c):
            return c[0] if isinstance(c, tuple) and len(c) == 1 else c
        col_keys = [_key(c) for c in cols]
        sampled_set = set(sampled)
        keep_cols = [c for c, k in zip(cols, col_keys) if k in sampled_set]
        dropped = [c for c, k in zip(cols, col_keys) if k not in sampled_set]
        if not dropped:
            return u_nk, []
        keep_keys = set(_key(c) for c in keep_cols)
        u_nk = u_nk[keep_cols]
        mask = u_nk.index.get_level_values("lambdas").isin(keep_keys)
        u_nk = u_nk.loc[mask]
        return u_nk, dropped

    @staticmethod
    def _sanitize_dHdl(dHdl):
        """Drop blown-up dV/dl frames before TI (final39).

        The MBAR u_nk filter never protected the TI headline: the dV/dl stream
        went straight into TI(). At a fully-decoupled vdw endpoint (lambda=1.000)
        pmemd can emit an endpoint-singularity dV/dl ('Energy at 0.0000 = ****'
        -> non-finite, or values of hundreds-thousands of kcal/mol). A single
        such frame poisons the TI integral (observed estimator spread ~1800
        kcal/mol -> dG_bind +95). We remove, PER LAMBDA:
          (a) non-finite frames, and
          (b) robust outliers: |x - median| > _DVDL_OUTLIER_MAD * MAD, also
              capped at an absolute |dV/dl| <= _DVDL_SANE_MAX_KCAL.
        Returns (clean_dHdl, n_dropped). Never raises (best-effort)."""
        try:
            import numpy as np, pandas as pd
        except ImportError:
            return dHdl, 0
        if dHdl is None or len(dHdl) == 0:
            return dHdl, 0
        col = dHdl.columns[0]
        x = dHdl[col].to_numpy(dtype=float, copy=True)
        finite = np.isfinite(x)
        keep = finite & (np.abs(np.where(finite, x, 0.0)) <= _DVDL_SANE_MAX_KCAL)
        # per-lambda robust outlier removal (group by all index levels except time)
        try:
            lvls = [l for l in dHdl.index.names if l != "time"]
            if lvls:
                grp = dHdl.reset_index().groupby(lvls)[col]
                med = grp.transform("median").to_numpy()
                mad = grp.transform(lambda s: (s - s.median()).abs().median()).to_numpy()
                mad = np.where(mad <= 0, np.nan, mad)
                robust = np.abs(x - med) <= (_DVDL_OUTLIER_MAD * 1.4826 * mad)
                robust = np.where(np.isnan(robust), True, robust).astype(bool)
                keep = keep & robust
        except Exception:  # noqa: BLE001
            pass
        n = int((~keep).sum())
        return (dHdl if n == 0 else dHdl.loc[keep]), n

    @staticmethod
    def _block_error(estimator_cls, data, kT, nblocks=_ERROR_BLOCKS):
        try:
            import numpy as np
            n = len(data)
            if n < nblocks * 5:
                return None
            edges = [int(round(k*n/nblocks)) for k in range(nblocks+1)]
            vals = []
            for k in range(nblocks):
                chunk = data.iloc[edges[k]:edges[k+1]]
                if len(chunk) < 3:
                    continue
                est = estimator_cls(); est.fit(chunk)
                vals.append(float(est.delta_f_.iloc[0,-1]*kT))
            if len(vals) < 2:
                return None
            return float(np.std(vals, ddof=1))
        except Exception:
            return None

    @staticmethod
    def _overlap_summary(mbar_estimator):
        """v2.4.25/26 (C): min nearest-neighbour overlap from a fitted
        alchemlyb MBAR (alchemlyb>=2.1 exposes .overlap_matrix; pymbar4
        backend). v2.4.26: log the REAL reason on failure instead of
        silently returning None, so 'overlap: null' is explainable."""
        import numpy as np
        O = None
        try:
            O = mbar_estimator.overlap_matrix
        except Exception as e1:
            try:
                mb = getattr(mbar_estimator, '_mbar', None)
                O = mb.compute_overlap()['matrix'] if mb is not None else None
            except Exception as e2:
                log.warning("  overlap matrix unavailable (%s / %s)", e1, e2)
                return None
        if O is None:
            log.warning("  overlap matrix is None (MBAR may not have converged)")
            return None
        try:
            O = np.asarray(O, dtype=float); k = O.shape[0]
            nn = [float(O[j, j+1]) for j in range(k-1)]
            weak = [j for j, v in enumerate(nn) if v < 0.03]
            return {"min_neighbor_overlap": (min(nn) if nn else None),
                    "n_weak_adjacent_pairs": len(weak),
                    "weak_pair_indices": weak}
        except Exception as e:
            log.warning("  overlap summary failed: %s", e)
            return None

    @staticmethod
    def _gauss_quad_ti(dvdl_by_lambda, nodes, weights):
        try:
            vals = []
            for nd in nodes:
                key = min(dvdl_by_lambda, key=lambda L: abs(L - nd))
                if abs(key - nd) > 1e-3:
                    return None
                vals.append(dvdl_by_lambda[key])
            return float(sum(w*v for w, v in zip(weights, vals)))
        except Exception:
            return None

    def run(self):
        try:
            import pandas as pd
            from alchemlyb.parsing.amber import extract_u_nk, extract_dHdl
            from alchemlyb.estimators import TI, BAR, MBAR
            from alchemlyb.preprocessing import statistical_inefficiency
        except ImportError as e:
            log.warning("alchemlyb not available (%s); trapezoid TI.", e)
            return self._finalize(self._fallback_ti())

        u_nk_list, dHdl_list, found = [], [], []
        _short_windows = []
        for lam in self.lambdas:
            out = self.dir / f"lambda_{lam:.3f}" / "prod.out"
            if not out.exists():
                log.warning("  missing %s", out); continue
            if not self._window_complete(self.dir, lam):
                _short_windows.append(round(float(lam), 3))
                log.warning("  lambda=%.3f prod.out present but INCOMPLETE "
                            "(no pmemd completion marker) -- excluded from "
                            "analysis as a missing window", lam)
                continue
            try:
                dHdl_list.append(extract_dHdl(out, T=self.T))
                u_nk_list.append(extract_u_nk(out, T=self.T))
                found.append(lam)
            except Exception as e:
                log.warning("  alchemlyb parse failed for lambda=%.3f: %s", lam, e)

        if not dHdl_list: return self._finalize(self._fallback_ti())

        import pandas as pd
        # v2.5.14 FIX 1: RAGGED MBAR-grid guard. Each window's prod.out records
        # MBAR energies on whatever lambda schedule was active when THAT window
        # ran. If a leg mixes schedules (e.g. 23 pre-2.5.11 windows + 5 refined
        # 28-grid windows), the per-window u_nk frames have DIFFERENT column sets;
        # pd.concat then yields a ragged/non-square matrix that MBAR/BAR cannot
        # solve (DLASCL / "lambda state not represented"). UNFIXABLE post-hoc --
        # the missing cross-energies were never computed -- so we fail loudly with
        # an actionable message instead of emitting silent garbage.
        _grids = [tuple(u.columns) for u in u_nk_list]
        if len(set(_grids)) > 1:
            from collections import Counter
            _sizes = Counter(len(g) for g in _grids)
            log.error("  INCONSISTENT MBAR grid across windows in %s: windows "
                      "report %d different state-count(s) %s. This leg mixes "
                      "lambda schedules and CANNOT be combined by MBAR/BAR. "
                      "Re-run the ENTIRE leg on a single schedule (config "
                      "vdw_lambdas/decharge_lambdas). TI (per-window dV/dl) is "
                      "unaffected and is reported below.",
                      self.dir, len(set(_grids)), dict(_sizes))
            return self._finalize(self._fallback_ti())
        dHdl = pd.concat(dHdl_list); u_nk = pd.concat(u_nk_list)
        u_nk, _n_masked = self._sanitize_u_nk(u_nk)
        dHdl, _n_dvdl_masked = self._sanitize_dHdl(dHdl)
        if _n_dvdl_masked:
            log.warning("  TI: dropped %d blown-up dV/dl frame(s) "
                        "(endpoint singularity guard)", _n_dvdl_masked)
        if _n_masked:
            log.warning("  sanitized u_nk: clipped %d soft-core end-state "
                        "singularity cell(s) to +%g kT (v2.5.14 clip guard; "
                        "zero Boltzmann weight, no rows/cols dropped)",
                        _n_masked, _UNK_SANE_MAX_KT)

        # final40 FIX 1: the alchemlyb amber parser restarts each window's
        # `time` index at 0, so pd.concat across lambdas yields DUPLICATE
        # (time) values. statistical_inefficiency then raises
        # "Duplicate time values found ... single, contiguous, sorted
        # timeseries", the decorrelation is skipped, and MBAR/BAR receive
        # malformed u_nk -> SVD divergence -> BAR ~ -1500 (the false 'spread').
        # Fix: decorrelate PER lambda group (dedup + sort each window), then
        # concat. This is the alchemlyb-recommended subsampling granularity.
        def _decorrelate_per_lambda(df):
            try:
                groups = []
                for _lam, g in df.groupby(level="lambdas"):
                    g = g[~g.index.duplicated(keep="first")].sort_index()
                    try:
                        g = statistical_inefficiency(g, g.iloc[:, 0])
                    except Exception as ee:  # noqa: BLE001
                        log.warning("  statistical_inefficiency skipped for "
                                    "lambda=%s (%s); using de-duplicated frames",
                                    _lam, ee)
                    groups.append(g)
                return pd.concat(groups) if groups else df
            except Exception as e:  # noqa: BLE001
                log.warning("  per-lambda decorrelation failed (%s); "
                            "using raw frames", e)
                return df
        dHdl = _decorrelate_per_lambda(dHdl)
        u_nk = _decorrelate_per_lambda(u_nk)

        kT = 0.001987204 * self.T
        estimators = {}
        _failed_estimators = {}   # final40 FIX 2: name -> reason (raised or diverged)
        overlap_info = None
        # final60 BUG 2 FIX (guard): if a whole lambda window contributed no
        # surviving frames, the sampled-lambda index no longer matches the
        # mbar_lambda state columns. Reconcile to the intersection so BAR/MBAR
        # receive a SQUARE reduced-potential matrix instead of raising
        # "Shape of passed values is (n,n), indices imply (n+1,n+1)".
        u_nk, _dropped_states = self._reconcile_u_nk_states(u_nk)
        if _dropped_states:
            log.warning("  u_nk reconciled: %d lambda-state column(s) had no "
                        "surviving samples and were excluded from MBAR/BAR "
                        "(%s). TI is unaffected.",
                        len(_dropped_states), _dropped_states)
        # v2.5.15 RANK GUARD: validate u_nk supports MBAR before solving.
        _mbar_rank, _mbar_cond, _mbar_nstates, _mbar_well_posed = \
            self._u_nk_rank_cond(u_nk)
        if not _mbar_well_posed:
            _failed_estimators["MBAR"] = (
                "ill-posed: u_nk rank %s of %s (cond=%.3e) -- under-determined "
                "MBAR; using BAR/TI. Densify vdw (use_dense_vdw=True) or sample "
                "longer." % (_mbar_rank, _mbar_nstates, _mbar_cond))
            log.warning("  MBAR ill-posed: u_nk rank %s of %s (cond=%.3e). "
                        "MBAR NOT computed; BAR/TI stand.",
                        _mbar_rank, _mbar_nstates, _mbar_cond)
        _cls = {"TI": TI, "BAR": BAR, "MBAR": MBAR}
        _plan = (("TI",TI(),dHdl),("BAR",BAR(),u_nk))
        if _mbar_well_posed:
            _plan = _plan + (("MBAR",MBAR(),u_nk),)
        for name, est, data in _plan:
            try:
                est.fit(data)
                rec = {
                    "dG_kcal_mol":  float(est.delta_f_.iloc[0,-1]*kT),
                    "err_kcal_mol": float(est.d_delta_f_.iloc[0,-1]*kT),
                }
                be = self._block_error(_cls[name], data, kT)
                if be is not None:
                    rec["block_err_kcal_mol"] = be
                estimators[name] = rec
                if name == "MBAR":
                    overlap_info = self._overlap_summary(est)
            except Exception as e:
                _failed_estimators[name] = "fit raised: %s" % e
                log.warning("  %s failed: %s", name, e)
        if overlap_info is not None:
            mno = overlap_info.get("min_neighbor_overlap")
            if mno is not None and mno < 0.03:
                log.warning("  LOW MBAR overlap: min nearest-neighbour %.4f "
                            "(%d weak adjacent pair(s)) -- add windows / sample longer.",
                            mno, overlap_info.get("n_weak_adjacent_pairs", 0))

        comp = self._completeness(found)
        if not comp["complete"]:
            log.warning("  INCOMPLETE leg %s: %d/%d windows analyzed; "
                        "missing %s. dG is reported but NOT trusted as a "
                        "headline result.", self.dir.name, comp["n_windows"],
                        comp["n_requested"], comp["missing_windows"])

        # --- v2.4.25 (B, OPTIONAL): Gaussian-quadrature TI override ---
        if getattr(self, "_gq_nodes", None) and getattr(self, "_gq_weights", None):
            gq = self._gauss_quad_ti(self._collect_dvdl(), self._gq_nodes, self._gq_weights)
            if gq is not None:
                estimators["TI"] = {"dG_kcal_mol": gq, "err_kcal_mol": float("nan"),
                                    "method": "gaussian_quadrature"}
                log.info("  TI via Gaussian quadrature (%d nodes): %+0.3f", len(self._gq_nodes), gq)
            else:
                log.warning("  gaussian_quadrature requested but dV/dl not at nodes; using standard TI.")

        # --- v2.4.21: TI-anchored headline selection ---
        # TI is the most robust estimator at decoupled endpoints (BAR/MBAR
        # overlap matrices are fragile). Reject BAR/MBAR if they disagree
        # with TI by more than _ESTIMATOR_CONSISTENCY_KCAL.
        ti = estimators.get("TI")
        ti_dG = ti.get("dG_kcal_mol") if ti else None

        def _consistent_with_ti(cand):
            if not self._is_sane(cand):
                return False
            if ti_dG is None:
                return True
            try:
                return abs(cand["dG_kcal_mol"] - ti_dG) <= _ESTIMATOR_CONSISTENCY_KCAL
            except (TypeError, ValueError):
                return False

        headline_name, headline = None, None
        for name in ("MBAR", "BAR", "TI"):
            cand = estimators.get(name)
            if name == "TI":
                if self._is_sane(cand):
                    headline, headline_name = cand, "TI"
                break
            if _consistent_with_ti(cand):
                headline, headline_name = cand, name
                break
        if headline is None:
            for name in ("TI", "BAR", "MBAR"):
                cand = estimators.get(name)
                if self._is_sane(cand):
                    headline, headline_name = cand, name
                    break
        if headline is None:
            return self._finalize(self._fallback_ti())

        for nm in ("MBAR", "BAR"):
            e = estimators.get(nm)
            if e is not None and not _consistent_with_ti(e):
                log.warning("  %s rejected: dG=%+.2f disagrees with TI=%+.2f "
                            "(tol %.1f) -> using %s", nm,
                            (e.get("dG_kcal_mol") or float('nan')),
                            (ti_dG if ti_dG is not None else float('nan')),
                            _ESTIMATOR_CONSISTENCY_KCAL, headline_name)

        # final40 FIX 2: an estimator that RAISED (SVD/DLASCL) or that
        # DIVERGED (|dG - TI| absurdly large, e.g. BAR=-1551 vs TI=-46 from a
        # collapsed u_nk) must NOT define the spread or condemn the leg as
        # non-converged. We exclude such estimators and record them explicitly.
        _DIVERGENCE_KCAL = 100.0   # |est - TI| beyond this == numerical failure
        spreads = {}
        for nm in ("BAR", "MBAR"):
            if nm not in estimators or ti_dG is None:
                continue
            val = estimators[nm].get("dG_kcal_mol")
            if val is None:
                continue
            d = val - ti_dG
            if abs(d) > _DIVERGENCE_KCAL:
                _failed_estimators[nm] = (
                    "diverged: |%s - TI| = %.1f kcal/mol (> %.0f) -- "
                    "collapsed u_nk / SVD failure; excluded from spread"
                    % (nm, abs(d), _DIVERGENCE_KCAL))
                estimators[nm]["diverged"] = True
                log.warning("  %s DIVERGED (dG=%.1f vs TI=%.1f); excluded "
                            "from spread/convergence.", nm, val, ti_dG)
                continue
            spreads[nm] = d
        max_spread = max((abs(v) for v in spreads.values()), default=0.0)
        # If every cross-energy estimator failed/diverged, the spread is
        # undefined -- TI stands alone. Report convergence on TI's own merits
        # (completeness handled separately) rather than a phantom huge spread.
        _only_ti = (len(spreads) == 0)

        # v2.5.36 RELIABILITY GATE
        _DVDL_SANE_MAX = 200.0
        _BARTI_SANE_MAX = 50.0
        _dvdl_now = self._collect_dvdl()
        _dvdl_max = max((abs(v) for v in _dvdl_now.values()), default=0.0)
        _reliable = True; _reasons = []
        if _dvdl_max > _DVDL_SANE_MAX:
            _reliable = False; _reasons.append(f"max|dV/dl|={_dvdl_max:.0f}>{_DVDL_SANE_MAX:.0f}")
        _bar = estimators.get("BAR"); _ti2 = estimators.get("TI")
        if _bar and _ti2:
            try:
                _d = abs(_bar["dG_kcal_mol"] - _ti2["dG_kcal_mol"])
                if _d > _BARTI_SANE_MAX:
                    _reliable = False; _reasons.append(f"|BAR-TI|={_d:.0f}>{_BARTI_SANE_MAX:.0f}")
            except (TypeError, ValueError, KeyError):
                pass
        if not _reliable:
            log.error("  RELIABILITY GATE FAILED for %s: %s. dG UNRELIABLE.",
                      self.dir.name, "; ".join(_reasons))
        # An incomplete leg must not advertise a trusted headline estimator.
        trusted_headline = headline_name if (comp["complete"] and _reliable) else None
        _leg_failed = ((not comp["complete"]) and (not _mbar_well_posed)
                       and ("BAR" in _failed_estimators))
        if _leg_failed:
            log.error("  LEG FAILED: %s incomplete (%d/%d; missing %s), MBAR "
                      "ill-posed and BAR diverged. dG is NOT reliable.",
                      self.dir.name, comp["n_windows"], comp["n_requested"],
                      comp["missing_windows"])
        return self._finalize({
            "dG_kcal_mol":           headline["dG_kcal_mol"],
            "estimators":            estimators,
            "dvdl":                  self._collect_dvdl(),
            "temperature":           self.T,
            "headline_estimator":    trusted_headline,
            "dG_reliable":            _reliable,
            "reliability_reasons":    _reasons,
            "max_abs_dvdl_kcal":      _dvdl_max,
            "estimator_used":        headline_name,
            "estimator_spread_kcal": max_spread,
            "estimator_converged":   max_spread <= _ESTIMATOR_CONSISTENCY_KCAL,
            "failed_estimators":     _failed_estimators,
            "only_ti_survived":      _only_ti,
            "mbar_u_nk_rank":        _mbar_rank,
            "mbar_u_nk_cond":        _mbar_cond,
            "mbar_n_states":         _mbar_nstates,
            "mbar_well_posed":       _mbar_well_posed,
            "overlap":               overlap_info,
            "block_err_kcal_mol":    headline.get("block_err_kcal_mol"),
            "short_windows":         _short_windows,
            "leg_failed":            _leg_failed,
            "dG_reliable":           bool(comp["complete"] and not _leg_failed),
            **comp,
        })

    def _finalize(self, result):
        corr_file = self.dir / "boresch_correction.txt"
        if corr_file.exists() and result.get("dG_kcal_mol") is not None:
            try:
                corr = float(corr_file.read_text().strip())
                result["dG_boresch_correction"] = corr  # v2.4.24 applied once in _leg()
                log.info("  Boresch correction recorded (applied once in _leg): %+0.3f", corr)
            except Exception as e:
                log.warning("  could not parse %s: %s", corr_file, e)
        complete = result.get("complete", True)
        missing  = result.get("missing_windows", [])
        incomplete_note = ("" if complete else
            f"# INCOMPLETE: {result.get('n_windows','?')}/"
            f"{result.get('n_requested','?')} windows; missing {missing}; "
            f"results below are UNRELIABLE\n")
        leg_label = self.dir.name   # v2.4.23: e.g. complex_decharge
        if result.get("dvdl"):
            d = result["dvdl"]
            (self.dir / "dvdl_summary.csv").write_text(
                f"# leg: {leg_label}\n" +
                incomplete_note +
                "leg,lambda,dvdl_kcal_mol\n" +
                "\n".join(f"{leg_label},{k:.3f},{v:.6f}"
                           for k, v in sorted(d.items())) + "\n")
        if result.get("estimators"):
            import csv
            head = result.get("headline_estimator")   # None if incomplete
            with open(self.dir / "dG_estimators.csv", "w", newline="") as f:
                f.write(f"# leg: {leg_label}\n")
                if incomplete_note:
                    f.write(incomplete_note)
                w = csv.writer(f)
                w.writerow(["leg","estimator","dG_kcal_mol","err_kcal_mol",
                            "headline","complete"])
                for name, v in result["estimators"].items():
                    err = v["err_kcal_mol"]
                    w.writerow([leg_label, name, f"{v['dG_kcal_mol']:.4f}",
                                "" if err is None else f"{err:.4f}",
                                "*" if (name == head and complete) else "",
                                "yes" if complete else "no"])
        if result.get("dvdl"):
            try:
                import matplotlib; matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                d = result["dvdl"]; xs = sorted(d); ys = [d[x] for x in xs]
                plt.figure(); plt.plot(xs, ys, "o-")
                plt.xlabel("lambda"); plt.ylabel("<dV/dlambda> (kcal/mol)")
                dG = result.get("dG_kcal_mol")
                title = f"{self.dir.name}"
                if dG is not None:
                    title += (f"   dG={dG:+.2f} kcal/mol "
                              f"[{result.get('headline_estimator', '?')}]")
                plt.title(title); plt.grid(alpha=0.3); plt.tight_layout()
                plt.savefig(self.dir / "dvdl.png", dpi=140); plt.close()
            except ImportError: pass
        return result


def relative_binding_dG(c, s):
    a = c.get("dG_kcal_mol"); b = s.get("dG_kcal_mol")
    return None if a is None or b is None else a - b


def absolute_binding_dG_from_legs(absolute_leg, solvent_leg):
    dG_c = absolute_leg.get("dG_kcal_mol")
    dG_s = solvent_leg.get("dG_kcal_mol")
    if dG_c is None or dG_s is None: return None
    return -(dG_c - dG_s)


def absolute_binding_dG(absolute_leg):
    log.warning("absolute_binding_dG(leg) is deprecated; "
                "use absolute_binding_dG_from_legs(absolute, solvent).")
    return absolute_leg.get("dG_kcal_mol")