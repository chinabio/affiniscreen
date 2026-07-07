"""Site / environment configuration -- the single source of truth for the
things that change when you move this workflow to a new machine or cluster.

Resolution order (first found wins):
    1. $BFEP_SITE_CONFIG                (explicit path, file)
    2. ~/.bfep/site_config.yaml|json    (per-user override)
    3. <repo>/site_config.yaml|json     (checked-in site default)
    4. built-in dataclass defaults      (this file)

The on-disk format may be YAML or JSON. YAML is used for writing when PyYAML is
installed; otherwise JSON is written (no hard dependency on PyYAML). A JSON file
is always valid YAML, so the loader accepts both regardless.

Nothing here imports Streamlit, so it is safe to use from the CLI / engine code.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict, fields
from pathlib import Path
from typing import List, Optional

ENV_VAR = "BFEP_SITE_CONFIG"
USER_DIR = Path.home() / ".bfep"
_REPO_ROOT = Path(__file__).resolve().parent.parent
_BASENAMES = ("site_config.yaml", "site_config.yml", "site_config.json")


# ---------------------------------------------------------------------------
# Dataclasses (defaults mirror the historical the cluster/the login node HPCConfig values)
# ---------------------------------------------------------------------------
@dataclass
class SchedulerConfig:
    # "lsf" is the only fully-supported backend today. "slurm" and "local" are
    # recognized so the GUI/scheduler layer can present them, but "slurm" is a
    # stub (see amber_md.schedulers) until validated on a real SLURM cluster.
    type: str = "lsf"
    # Command names (override if your site wraps them, e.g. "bsub" -> "lsf_bsub").
    submit_cmd: str = "bsub"
    query_cmd: str = "bjobs"
    cancel_cmd: str = "bkill"
    # Defaults surfaced as GUI field defaults.
    gpu_queue: str = "gpu"
    cpu_queue: str = "normalQ"
    project: str = "your-project"   # LSF project / SLURM account
    walltime: str = "24:00"
    max_concurrent: int = 8


@dataclass
class ResourceConfig:
    n_gpu: int = 1
    n_gpu_slots: int = 1
    gpu_request_mode: str = "rusage"          # "rusage" | "slots"
    gpu_cpu_cores: int = 8                     # cores requested alongside 1 GPU
    fep_gpu_cores: int = 16                    # cores on a FEP window job
    fep_mem_mb: int = 8192
    cpu_settle_cap_gpu: int = 32               # AMBERMD_CPU_NP cap on GPU nodes
    cpu_settle_cap_cpu: int = 48               # ... on the CPU queue
    n_cpu: int = 8                             # MM/GBSA MPI ranks default
    avoid_hosts: List[str] = field(default_factory=list)


@dataclass
class AmberConfig:
    # One of these is enough. amberhome is used to locate binaries; module_load
    # is emitted into generated job scripts; sander/pmemd override explicit paths.
    amberhome: str = ""                         # $AMBERHOME (empty = rely on env)
    module_load: List[str] = field(
        default_factory=lambda: ["gcc/11.5", "cuda/11.8", "amber/22.8"])
    pmemd_cuda: str = "pmemd.cuda"              # binary name or absolute path
    charge_method: str = "bcc"


@dataclass
class OpenFEConfig:
    # How to reach the OpenFE/OpenMM env from the login/compute node.
    conda_env: str = "openfe"                   # `conda run -n <env>` if set
    python_bin: str = ""                        # explicit python (overrides env)
    small_molecule_ff: str = "openff-2.1.0"


@dataclass
class PathsConfig:
    default_run_dir: str = "~/Run_dir"
    scratch_dir: str = ""                        # empty = same as run dir
    venv_activate: str = "./activate_amber_md.sh"


@dataclass
class SiteConfig:
    scheduler: SchedulerConfig = field(default_factory=SchedulerConfig)
    resources: ResourceConfig = field(default_factory=ResourceConfig)
    amber: AmberConfig = field(default_factory=AmberConfig)
    openfe: OpenFEConfig = field(default_factory=OpenFEConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)

    # ---- (de)serialization -------------------------------------------------
    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SiteConfig":
        d = d or {}
        def _sub(kls, key):
            raw = d.get(key, {}) or {}
            valid = {f.name for f in fields(kls)}
            return kls(**{k: v for k, v in raw.items() if k in valid})
        return cls(
            scheduler=_sub(SchedulerConfig, "scheduler"),
            resources=_sub(ResourceConfig, "resources"),
            amber=_sub(AmberConfig, "amber"),
            openfe=_sub(OpenFEConfig, "openfe"),
            paths=_sub(PathsConfig, "paths"),
        )


# ---------------------------------------------------------------------------
# Load / save
# ---------------------------------------------------------------------------
def _parse(text: str) -> dict:
    """Parse YAML if available, else JSON. JSON is a subset of YAML, so JSON
    files load either way."""
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text) or {}
    except Exception:
        return json.loads(text or "{}")


def config_search_paths() -> list[Path]:
    paths: list[Path] = []
    env = os.environ.get(ENV_VAR)
    if env:
        paths.append(Path(env).expanduser())
    for base in _BASENAMES:
        paths.append(USER_DIR / base)
    for base in _BASENAMES:
        paths.append(_REPO_ROOT / base)
    return paths


def resolve_config_path() -> Optional[Path]:
    for p in config_search_paths():
        if p.is_file():
            return p
    return None


def load() -> SiteConfig:
    p = resolve_config_path()
    if p is None:
        return SiteConfig()
    try:
        return SiteConfig.from_dict(_parse(p.read_text()))
    except Exception:
        # Never let a malformed site file break the app; fall back to defaults.
        return SiteConfig()


def default_save_path() -> Path:
    """Where `save()` writes when no explicit path is given: per-user file."""
    env = os.environ.get(ENV_VAR)
    if env:
        return Path(env).expanduser()
    return USER_DIR / "site_config.yaml"


def save(cfg: SiteConfig, path: Optional[os.PathLike | str] = None) -> Path:
    out = Path(path).expanduser() if path else default_save_path()
    out.parent.mkdir(parents=True, exist_ok=True)
    data = cfg.to_dict()
    try:
        import yaml  # type: ignore
        text = yaml.safe_dump(data, sort_keys=False, default_flow_style=False)
        if out.suffix == ".json":
            text = json.dumps(data, indent=2)
    except Exception:
        text = json.dumps(data, indent=2)
        if out.suffix in (".yaml", ".yml"):
            out = out.with_suffix(".json")
    out.write_text(text)
    return out


# Convenience: a process-wide cached instance (re-read with load()).
_cached: Optional[SiteConfig] = None

def get() -> SiteConfig:
    global _cached
    if _cached is None:
        _cached = load()
    return _cached

def refresh() -> SiteConfig:
    global _cached
    _cached = load()
    return _cached
