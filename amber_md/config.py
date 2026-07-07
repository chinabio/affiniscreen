"""Centralized configuration (v2.4.20)."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional
import json


@dataclass
class SystemConfig:
    protein_ff: str = "ff19SB"
    ligand_ff: str = "gaff2"
    water_model: str = "tip3p"
    # v2.5.71: 16.0 -> 12.0. The 16 A buffer was set in v2.4.21 to stop a
    # 'GPU small-box crash on small ligands' -- but that instability was a
    # SOFT-CORE / box-drift problem, root-caused and fixed in v2.5.70
    # (stage-aware hard-vdW decharge) together with the 44-window vdw
    # schedule and the MD-only NaN guard. With those fixes the oversized box
    # is no longer needed. 12 A matches the validated OpenFE default and
    # leaves safe headroom over the 10 A nonbonded cutoff for a decoupling
    # (reorienting) soft-core ligand; 10 A would be uncomfortably tight.
    # On this 148k-atom test system (solute span ~93 A, 91% water) 12 A
    # cuts ~16 ns/day -> ~21 ns/day with no accuracy change. Validate a
    # single ligand before committing a full campaign.
    box_buffer_A: float = 12.0
    box_shape: str = "octahedral"
    neutralize: bool = True
    salt_conc_M: float = 0.15
    ligand_charge: int = 0
    ligand_multiplicity: int = 1
    charge_method: str = "bcc"
    ion_method: str = "rand"


@dataclass
class MDConfig:
    min_maxcyc: int = 10000
    min_ncyc: int = 5000
    min_restraint_wt: float = 5.0
    heat_nsteps: int = 50000
    heat_T_start: float = 100.0
    heat_T_end: float = 298.0
    equil_nsteps: int = 500000
    equil_restraint_wt: float = 1.0
    posres_mask: str = ""   # v2.5.19: solute heavy-atom mask held through
                            # dens+eq (e.g. "!:WAT,Na+,Cl-,K+ & !@H="); 
                            # "" keeps the prior ntr=0 equilibration.
    prod_nsteps: int = 10_000_000   # final52: 50 ns -> 10 ns (unified default)
    prod_dt_ps: float = 0.001
    prod_print_freq: int = 5000
    temperature_K: float = 298.0
    pressure_bar: float = 1.0
    cutoff_A: float = 10.0
    ntc: int = 2
    ntf: int = 2


@dataclass
class HPCConfig:
    """the login node-specific defaults (the cluster GPU nodes: 64 cores / 4 GPUs each).

    v2.4.5 changes:
      * n_cpu default = 8. With 64 cores and 4 GPUs per the GPU queue node, the
        per-GPU fair share is 16 cores; 8 leaves headroom for the OS and for
        a second GPU job to co-tenant the node. Drives mpirun -np in the
        in-job MM/GBSA step (which actually reads $LSB_DJOB_NUMPROC at
        runtime so the rank count always matches LSF's grant).
      * mmgbsa_n_cpu: cap on MM/GBSA MPI ranks, decoupled from job slot
        count. None means "use n_cpu". Lets you e.g. request 8 slots but
        only use 4 for MM/GBSA if you want lower memory pressure.
    """
    queue_gpu: str = "gpu"
    queue_cpu: str = "normalQ"
    project: str = "your-project"
    job_name: str = "amberMD"
    walltime: str = "24:00"
    n_cpu: int = 8
    n_gpu: int = 1
    # final53: LSF slots requested on the GPU queue. On the cluster `gpu`, slots
    # map to GPUs, so this MUST default to n_gpu (1) -- not n_cpu -- or a
    # single-GPU pmemd.cuda job grabs 8 GPUs (8*the GPU queue) and idles 7.
    n_gpu_slots: int = 1
    # final54: request extra CPU cores ALONGSIDE the 1 GPU for parallel
    # MM/GBSA scoring (MMPBSA.py.MPI). gpu_cpu_cores = cores to request via
    # `-n` while the GPU is requested as a RESOURCE (rusage[ngpus_physical])
    # so cores do NOT inflate the GPU count on queues where slots==GPUs.
    #   gpu_request_mode="slots"   -> legacy: -n=n_gpu_slots, -gpu num=N
    #                                 (use when slots==cores)
    #   gpu_request_mode="rusage"  -> -n=gpu_cpu_cores, GPU via
    #                                 rusage[ngpus_physical=n_gpu]
    #                                 (use when slots==GPUs, e.g. the cluster `gpu`)
    gpu_request_mode: str = "rusage"
    gpu_cpu_cores: int = 8
    mmgbsa_n_cpu: Optional[int] = None
    fep_mem_mb: int = 8192  # host RAM (MB) reserved per FEP window (the cluster mem enforcement)
    # v2.5.34: CPU cores requested (-n) on the GPU FEP array job, alongside
    # the 1 GPU. GPU MD uses ~1-2; the surplus fuels the rare inline CPU
    # density-settle fallback (mpirun -np). 16 is trivial on 32-64 core
    # the cluster GPU nodes. The settle further caps to 32 (GPU node) / 48 (normalQ).
    fep_gpu_cores: int = 16
    # v2.5.8: persistent GPU host blocklist. Hosts listed here are excluded
    # from EVERY generated GPU job via -R "select[hname!=...]" (initial launch
    # AND self-heal reruns). Set to nodes with known driver/ECC issues, e.g.
    # ("gpu-node-01", "gpu-node-02"). Empty tuple = no exclusion.
    avoid_hosts: tuple = ()
    modules: tuple = ("gcc/11.5", "cuda/11.8", "amber/22.8")
    # NOTE (v2.5.3): submit.py._header() resolves this to an ABSOLUTE
    # path at script-generation time (covers GUI/batch/run_amber.py).
    # The relative default is retained only for backward compatibility.
    venv_activate: Optional[str] = "./activate_amber_md.sh"


@dataclass
class MMGBSAConfig:
    enabled: bool = True
    igb: int = 8
    salt_conc: float = 0.15
    start_frame: int = 1
    end_frame: int = 0
    stride: int = 1
    decomposition: bool = False
    decomp_residues: str = ""


@dataclass
class FEPConfig:
    lambdas: tuple = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                      0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
    timask1: str = ":LIA"
    timask2: str = ":LIB"
    scmask1: str = ":LIA"
    scmask2: str = ":LIB"
    nstlim_eq: int = 250_000
    nstlim_prod: int = 20_000_000  # v2.5.75: 10->20 ns (10-14 day budget; double
    #                              # (decharge not converged at 2 ns); 5 ns/window
    # v2.5.64: the Option-A restraint leg is plain NPT equilibration for an
    # ANALYTIC Boresch correction (not sampling-limited), so it needs only
    # enough time to relax the restrained ensemble. 2 ns @1fs is a safety
    # margin over the ~1 ns that suffices; TI legs keep the full nstlim_prod.
    # Set to 0 to fall back to nstlim_prod.
    restraint_nstlim_prod: int = 2_000_000
    dt_ps: float = 0.001
    temperature_K: float = 298.0
    cutoff_A: float = 10.0
    ntpr: int = 1000
    ntwx: int = 5000
    # v2.5.16: AMBER-recommended soft pair (was 0.5/12.0, which overlapped
    # soft-core cores at mid-lambda and drove the 0.6-0.9 blow-ups).
    scalpha: float = 0.2
    scbeta: float = 50.0
    two_stage: bool = True
    crgmask: str = ":LIG"
    # v2.5.23 (Option B): dual-copy restraint-leg masks. The complex_restraint
    # leg runs on a dedicated  combine { LIG LIG protein }  topology where the
    # ligand is :1 (real) / :2 (dummy). Ordinary TI, ifsc=0, fixed-k Boresch.
    restraint_timask1: str = ":1"
    restraint_timask2: str = ":2"
    restraint_crgmask: str = ""            # v2.5.24: empty -> no crgmask (match FEP-SPell-ABFE)
    # build the dedicated dual-copy restraint topology before the leg (else the
    # leg would reuse complex.prmtop and hit the v2.5.20-2.5.22 failures).
    build_restraint_topology: bool = False   # v2.5.36 (Option A): single-copy lambda-scaled Boresch on real complex.prmtop
    restraint_reion: bool = False          # v2.5.24: inherit neutral equilibrated system
    # v2.4.25 (A): GTI softcore controls (BAT.py-proven).
    gti_add_sc: int = 1
    # v2.5.16: complete GTI soft-core control set (FEP-SPell-proven). The
    # ele/vdw/lam_sch trio was missing -> half-configured TI soft-core ->
    # stiff dV/dl -> mid-lambda divergence. gti_chg_keep 0 -> 1.
    gti_chg_keep: int = 1
    gti_lam_sch: int = 1
    gti_ele_sc: int = 1
    gti_vdw_sc: int = 1
    # v2.5.70: STAGE-AWARE decharge soft-core. ROOT CAUSE of the complex_decharge
    # analysis failure (TI=-56 vs BAR=-19, MBAR SVD non-convergent, **** / 4.8e8
    # entries in the high-lambda MBAR table): the decharge stage put :LIG in the
    # vdW soft-core region (scmask1=:LIG, gti_vdw_sc=1) AND kept soft-core charges
    # (gti_chg_keep=1) WHILE clambda removed those same charges with no crgmask.
    # Opposite partial charges then approached through a SOFTENED vdW core ->
    # -q.q/r charging singularity -> dV/dl ~ -507 kT and an unsolvable u_nk.
    # Fix: on the DECHARGE stage only, decharge with HARD vdW -- do not soften
    # vdW (gti_vdw_sc=0) and let charges be genuinely removed (gti_chg_keep=0).
    # The VDW stage is unchanged (charges already held off via crgmask=:LIG).
    decharge_gti_vdw_sc: int = 0      # decharge: HARD vdW (no soft-core on vdW)
    decharge_gti_chg_keep: int = 0    # decharge: actually remove charges
    gti_scale_beta: int = 0
    gti_cut_sc: int = 0
    gti_cut: int = 1
    tishake: int = 1
    logdvdl: int = 0
    # v2.5.16: dedicated Boresch restraint-removal leg schedule (16 windows,
    # published FEP-SPell grid). Used by the new complex_restraint leg.
    restraint_lambdas: tuple = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.8, 0.85,
                                0.9, 0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0)
    # v2.5.51: FRONT-LOADED restraint schedule (quadratic ramp). The default
    # FEP-SPell grid jumps 0.0 -> 0.15 -> 0.30, a single large step exactly where
    # the Boresch restraint switches on. lig_12944901 detonated at 0.15 with an
    # intermittent steric clash that NO restart strategy (v2.5.49 escalation incl.
    # CPU density-settle) could clear -- the signature of an UNDER-SAMPLED window.
    #
    # SYNTHESIZED from three independent sources, not one:
    #   * FEP+ (Schrodinger) production .msj for THIS EXACT ligand, which converged:
    #     its alchemical_boresch restraint-ON schedule is a smooth QUADRATIC ramp
    #     with a ~0.0009 first step and ~38% of windows below lambda=0.15.
    #   * BAT.py attach_rest (Heinzelmann & Gilson): front-load where restraint engages.
    #   * the reference platform/ByteDance zero-shot: smooth/uniform schedules, no per-system tuning;
    #     keep the cited default unchanged.
    # We reproduce the FEP+ SHAPE (quadratic, ~35% windows <0.15) at a TI-affordable
    # 20 windows (FEP+ used 34 cheap RE-umbrella steps). Spacing fixes the lambda=0.15
    # crash; per-window sim time is a separate lever.
    #
    # OPT-IN (use_fine_restraint_lambdas=True). Default restraint_lambdas above
    # remains the published FEP-SPell grid.
    restraint_lambdas_fine: tuple = (
        0.0, 0.004, 0.016, 0.036, 0.064, 0.1, 0.144, 0.196, 0.256, 0.324,
        0.4, 0.5, 0.6, 0.7, 0.8, 0.875, 0.925, 0.95, 0.975, 1.0,
    )
    use_fine_restraint_lambdas: bool = False
    # v2.4.26: gentle density-equilibration barostat (stops soft-core box
    # collapse at intermediate lambda). barostat=1 Berendsen, loose taup.
    dens_barostat: int = 1
    dens_taup: float = 5.0
    # final60 BUG 1: gentle EQUILIBRATION integrator. Equilibration used to
    # revert to dt=dt_ps (2 fs) + MC barostat off the gentle dens.rst, which
    # diverged on half-decoupled soft-core ligands at mid-lambda (T~15,000 K).
    # These knobs let _eq_in inherit the density-stage stability. Production is
    # unchanged (still MC barostat / dt_ps).
    eq_dt_ps: float = 0.001       # half step during equilibration
    eq_gamma_ln: float = 5.0      # stronger Langevin coupling
    eq_taup: float = 5.0          # loose barostat relaxation
    eq_barostat: int = 1          # 1 = Berendsen (stable) during eq
    eq_temp_max_K: float = 1000.0 # stability gate: fail eq above this peak T
    do_heat: bool = True
    heat_nstlim: int = 100_000
    heat_dt_ps: float = 0.001
    heat_gamma_ln: float = 5.0
    heat_T_start: float = 5.0
    heat_ramp_frac: float = 0.8
    heat_restraint_wt: float = 5.0
    heat_temp_max_K: float = 1000.0  # v2.5.26: stability gate for the heat stage
    vlimit: float = 20.0
    posres_mask_default: str = "!:WAT,HOH,Na+,Cl-,K+,Mg2+,Ca2+,Zn2+ & !@H="  # v2.5.28: no '*' wildcards
    # v2.4.25 (B, OPTIONAL, default OFF): Gaussian-quadrature TI.
    gaussian_quadrature: bool = False
    gq_nodes: tuple = (0.02544, 0.12923, 0.29707, 0.5,
                       0.70292, 0.87076, 0.97455)
    gq_weights: tuple = (0.06474, 0.13985, 0.19091, 0.20897,
                         0.19091, 0.13985, 0.06474)
    # v2.4.21: 9 -> 15 windows, denser at endpoints. The 9-window grid
    # left TI vs BAR ~17 kcal/mol apart (under-converged decharge dV/dl).
    # v2.5.70: 15 -> 21 windows. The old grid was too coarse across the
    # high-lambda electrostatic-decoupling singularity (dV/dl swung
    # +102 -> -507 -> -88 kcal/mol over 0.625->0.95). Even with the hard-vdW
    # decharge fix above, the residual charging curvature near full decoupling
    # needs finer spacing for stable TI trapezoid + adequate BAR/MBAR overlap.
    # Densified 0.5-0.975; endpoints/low-lambda unchanged. NOTE: changing this
    # schedule means any partially-completed decharge leg must be re-run from
    # scratch (MBAR needs one consistent lambda set). Pass --decharge-lambdas
    # to override (e.g. to reproduce the legacy 15-window grid).
    # v2.5.73: densified 0.725-0.925. The finished 21-window leg showed the
    # dV/dL well swing ~160 kcal/mol per window across 0.775-0.863 -> MBAR
    # unsolvable (6504 foreign-E > 0). Inserted 0.75, 0.80, 0.8438, 0.9125
    # to halve the steepest gaps (21 -> 26 windows).
    # v2.5.74: 26 -> 30 windows. With a week of 8-GPU budget we densify the
    # 0.72-0.93 charging well to <=0.0375 spacing everywhere (was up to 0.05),
    # to guarantee MBAR nearest-neighbour overlap across the steepest dV/dL.
    # v2.5.77: (1) rounded to 3 decimals -- prod.in writes mbar_lambda %.3f, so
    # 4-dp points truncated (0.5625->0.562 etc) and 11/30 windows produced no
    # u_nk -> MBAR dead. (2) DATA-DRIVEN inserts 0.656 + 0.706: the diagnostic
    # run showed dV/dl +49(lam0.5)->-93(lam0.688), slope ~-755/unit-lambda with
    # print overflows at the well bottom; halving the 0.625-0.725 gaps tames the
    # per-window perturbation and lifts overlap. 30 -> 32 windows.
    decharge_lambdas: tuple = (
        0.0, 0.025, 0.05, 0.1, 0.175, 0.275,
        0.375, 0.5, 0.562, 0.625, 0.656, 0.688,
        0.706, 0.725, 0.738, 0.75, 0.762, 0.775,
        0.788, 0.8, 0.812, 0.825, 0.844, 0.862,
        0.881, 0.9, 0.912, 0.925, 0.938, 0.95,
        0.975, 1.0)
    decharge_lambdas_legacy: tuple = (0.0, 0.025, 0.05, 0.1, 0.175,
                               0.275, 0.375, 0.5, 0.625, 0.725,
                               0.825, 0.9, 0.95, 0.975, 1.0)
    # final41: vdW decoupling schedule. The bare lambda=1.000 endpoint caused
    # an end-state softcore singularity in the smoke test. We KEEP full
    # decoupling (dropping 1.0 would bias dG by the missing 0.95->1.0 tail) but
    # add fine spacing 0.95/0.975/0.99/1.0 so each near-endpoint dV/dl step is
    # small and numerically stable. To run a hard 0.95 cap instead, pass an
    # explicit --vdw-lambdas (and accept the truncated-tail bias).
    # v2.5.11: refined vdW schedule. Under the original uniform 0.05 spacing,
    # complex_vdw windows at lambda 0.70/0.75 hit eq instability (exit 71) and
    # 0.15/.../0.70 hit prod box drift (exit 255) -- the classic soft-core
    # danger zone where the ~60-85% decoupled ligand core lets atoms nearly
    # overlap and dV/dl turns stiff. We halve the spacing across 0.6-0.85 so
    # each window's perturbation (and box response) is smaller and stable.
    # Endpoints/low-lambda are unchanged. NOTE: changing this schedule means
    # any partially-completed vdw leg must be re-run from scratch (MBAR needs a
    # single consistent lambda set). Pass --vdw-lambdas to override.
    # v2.5.11 28-window schedule. Kept as the default + fallback.
    vdw_lambdas_legacy: tuple = (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,
                          0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.625,
                          0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8,
                          0.825, 0.85, 0.9, 0.95, 0.975, 0.99, 1.0)
    # v2.5.16: published FEP-SPell vdW schedule (44 windows, ultra-dense
    # across 0.575-0.80 -- the real soft-core danger zone). This lifts
    # nearest-neighbour MBAR overlap above the 0.03 floor AND keeps each
    # window's perturbation small enough to stay numerically stable.
    vdw_lambdas: tuple = (
        0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.475,
        0.50, 0.525, 0.55, 0.575, 0.585, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64,
        0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75,
        0.76, 0.77, 0.78, 0.79, 0.80, 0.825, 0.85, 0.875, 0.90, 0.95, 1.0)
    # v2.5.15: denser 40-window schedule (danger zone 0.55-0.90 halved
    # to 0.025) to lift nearest-neighbour overlap / u_nk rank so MBAR
    # can be well-posed. Opt in with use_dense_vdw=True.
    # v2.5.74: 40 -> 50 windows. Halve spacing through the 0.55-0.90 soft-core
    # danger zone to <=0.0125 so u_nk rank stays full and MBAR is well-posed
    # (the decharge MBAR collapse warned this zone is overlap-critical for vdw too).
    # v2.5.77: rounded to 3 decimals (same %.3f mbar_lambda truncation bug).
    vdw_lambdas_dense: tuple = (
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25,
        0.3, 0.35, 0.4, 0.45, 0.475, 0.5,
        0.525, 0.55, 0.562, 0.575, 0.588, 0.6,
        0.612, 0.625, 0.638, 0.65, 0.662, 0.675,
        0.688, 0.7, 0.712, 0.725, 0.738, 0.75,
        0.762, 0.775, 0.788, 0.8, 0.812, 0.825,
        0.838, 0.85, 0.862, 0.875, 0.888, 0.9,
        0.912, 0.925, 0.938, 0.95, 0.962, 0.975,
        0.988, 1.0)
    use_dense_vdw: bool = True    # v2.5.73: default ON. vdw hits the same
    #                            # 0.55-0.90 overlap collapse as decharge did;
    #                            # dense 40-window grid halves spacing there.


@dataclass
class FEPWorkflowConfig:
    enabled: bool = False
    params: FEPConfig = field(default_factory=FEPConfig)
    complex_prmtop: Optional[Path] = None
    complex_inpcrd: Optional[Path] = None
    solvent_prmtop: Optional[Path] = None
    solvent_inpcrd: Optional[Path] = None


@dataclass
class WorkflowConfig:
    work_dir: Path = Path("./run")
    complex_pdb: Path = Path("complex.pdb")
    ligand_resname: str = "LIG"
    solvent_mask: str = ":WAT,K+,Na+,Cl-"   # v2.5.18 charge-correction
    water_model: str = "tip3p"               # v2.5.18 charge-correction
    ligand_input: Optional[Path] = None
    protein_input: Optional[Path] = None
    system: SystemConfig = field(default_factory=SystemConfig)
    md: MDConfig = field(default_factory=MDConfig)
    hpc: HPCConfig = field(default_factory=HPCConfig)
    mmgbsa: MMGBSAConfig = field(default_factory=MMGBSAConfig)
    fep: FEPWorkflowConfig = field(default_factory=FEPWorkflowConfig)
    submit: bool = True
    monitor: bool = True
    auto_protonation: bool = True
    protonation_overrides: Optional[dict] = None

    def save(self, path):
        def _enc(o):
            if isinstance(o, Path):
                return str(o)
            if isinstance(o, tuple):
                return list(o)
            return o
        Path(path).write_text(json.dumps(asdict(self), default=_enc, indent=2))

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text())
        fep_d = d.get("fep", {})
        params_d = fep_d.get("params", {})
        if "lambdas" in params_d:
            params_d = {**params_d, "lambdas": tuple(params_d["lambdas"])}
        params_d.setdefault("scalpha", 0.5)
        params_d.setdefault("scbeta", 12.0)
        params = FEPConfig(**params_d) if params_d else FEPConfig()
        fep = FEPWorkflowConfig(
            enabled=fep_d.get("enabled", False), params=params,
            complex_prmtop=Path(fep_d["complex_prmtop"]) if fep_d.get("complex_prmtop") else None,
            complex_inpcrd=Path(fep_d["complex_inpcrd"]) if fep_d.get("complex_inpcrd") else None,
            solvent_prmtop=Path(fep_d["solvent_prmtop"]) if fep_d.get("solvent_prmtop") else None,
            solvent_inpcrd=Path(fep_d["solvent_inpcrd"]) if fep_d.get("solvent_inpcrd") else None)

        # v2.4.5: tolerate older saved configs that lack mmgbsa_n_cpu.
        hpc_d = {**d["hpc"], "modules": tuple(d["hpc"]["modules"])}
        hpc_d.setdefault("mmgbsa_n_cpu", None)

        return cls(
            work_dir=Path(d["work_dir"]), complex_pdb=Path(d["complex_pdb"]),
            ligand_resname=d.get("ligand_resname", "LIG"),
            ligand_input=Path(d["ligand_input"]) if d.get("ligand_input") else None,
            protein_input=Path(d["protein_input"]) if d.get("protein_input") else None,
            system=SystemConfig(**d["system"]), md=MDConfig(**d["md"]),
            hpc=HPCConfig(**hpc_d),
            mmgbsa=MMGBSAConfig(**d["mmgbsa"]), fep=fep,
            submit=d.get("submit", True), monitor=d.get("monitor", True),
            # v2.5.0 FIX: these were silently dropped on load, so config-file
            # protonation settings were ignored. Now round-trip correctly.
            auto_protonation=d.get("auto_protonation", True),
            protonation_overrides=d.get("protonation_overrides", None))