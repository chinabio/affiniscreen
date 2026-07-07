"""Shared widgets and helpers for the multipage GUI."""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import os, re, subprocess
from pathlib import Path
from datetime import datetime
import streamlit as st

def list_dir(path, exts=None):
    try:
        p = Path(path).expanduser().resolve()
        if not p.is_dir(): return []
        out = []
        for f in sorted(p.iterdir()):
            if f.is_file():
                if exts is None or f.suffix.lower() in exts:
                    out.append(str(f))
        return out
    except Exception:
        return []

def file_picker(label, key, exts, default_dir="~"):
    col1, col2 = st.columns([3, 2])
    init = st.session_state.get(f"{key}_dir") or os.path.expanduser(default_dir)
    cur_dir = col1.text_input(f"{label} -- directory", init, key=f"{key}_dir")
    files = list_dir(cur_dir, exts)
    options = ["<pick a file>"] + [Path(f).name for f in files]
    picked = col2.selectbox(f"{label} -- file", options, key=f"{key}_pick")

    full_key = f"{key}_full"
    # If user picked a file from the dropdown, push that path into the override
    # field's session state BEFORE the widget is instantiated. Streamlit ignores
    # `value=` on widgets that already have a key, so we must use session_state.
    if picked != "<pick a file>":
        computed = str(Path(cur_dir).expanduser().resolve() / picked)
        existing = st.session_state.get(full_key, "")
        try:
            stale = (
                not existing
                or Path(existing).name != picked
                or str(Path(existing).parent) != str(Path(cur_dir).expanduser().resolve())
            )
        except Exception:
            stale = True
        if stale:
            st.session_state[full_key] = computed

    override = st.text_input(f"{label} -- full path (override)", key=full_key)
    return (override or "").strip() or None

def _subdirs(path):
    """Return sorted immediate sub-directories of *path* (silent on errors)."""
    try:
        p = Path(path).expanduser().resolve()
        if not p.is_dir():
            return []
        return sorted([d for d in p.iterdir() if d.is_dir()
                       and not d.name.startswith(".")],
                      key=lambda d: d.name.lower())
    except (PermissionError, OSError):
        return []


def dir_picker(label, key, default_dir="~", browse=True):
    """Directory chooser: type a path OR browse a folder tree.

    Backward-compatible: same signature/return (resolved path string or None).

    Streamlit rule respected: we NEVER assign to the text-input's own
    session_state key after the widget is created. Instead, browser buttons
    write to a separate ``<key>_pending`` slot and call st.rerun(); on the next
    run that pending value becomes the text_input's default *before* the widget
    is instantiated. This avoids
    'st.session_state.<key> cannot be modified after the widget ... instantiated'.
    """
    tkey = f"{key}_dir"           # text-input widget key (unchanged -> compatible)
    pend = f"{key}_pending"       # pending value to seed the text box next run
    bkey = f"{key}_browse_cwd"    # browser current directory

    # 1) Resolve the initial value BEFORE instantiating the widget.
    #    Streamlit forbids passing BOTH value= AND having the widget's own
    #    session_state key already set -- that emits:
    #      "The widget with key ... was created with a default value but also
    #       had its value set via the Session State API."
    #    So: if the key is already in session_state (set on a prior run or by a
    #    consumed _pending slot), let session_state be the SOLE source of truth
    #    and DO NOT pass value=. Only seed value= on the very first run.
    if pend in st.session_state:
        # A browser button selected a folder: seed the key before instantiation.
        st.session_state[tkey] = st.session_state.pop(pend)

    # 2) The text input (source of truth).
    if tkey in st.session_state:
        # Key already present -> session_state wins; no value= (avoids warning).
        val = st.text_input(label, key=tkey)
    else:
        # First run for this key -> provide the initial default via value=.
        init = os.path.expanduser(default_dir)
        val = st.text_input(label, value=init, key=tkey)
    p = Path(val).expanduser().resolve() if val else None

    # 3) Optional browser. Buttons only set the pending slot + rerun; they do
    #    NOT touch tkey directly (which is now illegal post-instantiation).
    if browse:
        with st.expander("Browse folders", expanded=False):
            if bkey not in st.session_state:
                st.session_state[bkey] = (str(p) if (p and p.is_dir())
                                          else os.path.expanduser("~"))
            cwd = Path(st.session_state[bkey]).expanduser().resolve()

            st.caption(f"Current: `{cwd}`")
            c0, c1, c2 = st.columns([1, 1, 2])
            if c0.button("Up", key=f"{key}_up", width="stretch"):
                st.session_state[bkey] = str(cwd.parent)
                st.rerun()
            if c1.button("Go to typed", key=f"{key}_sync", width="stretch"):
                if p and p.is_dir():
                    st.session_state[bkey] = str(p)
                    st.rerun()
            if c2.button("Select this folder", key=f"{key}_sel",
                         type="primary", width="stretch"):
                st.session_state[pend] = str(cwd)   # seed text box next run
                st.rerun()

            subs = _subdirs(cwd)
            if subs:
                names = [d.name for d in subs]
                pick = st.radio("Sub-folders", names, index=None,
                                key=f"{key}_subpick")
                if pick:
                    st.session_state[bkey] = str(cwd / pick)
                    st.rerun()
            else:
                st.caption("(no accessible sub-folders here)")

    # 4) Validate + return (value reflects the widget this run).
    if p and not p.exists():
        st.caption(f"Directory does not exist yet: `{p}`")
    elif p and not p.is_dir():
        st.caption(f"Not a directory: `{p}`")
    return str(p) if p else None

def get_lsf_jobs(user=None):
    user = user or os.environ.get("USER","")
    try:
        cp = subprocess.run(
            ["bjobs","-u",user,"-o","JOBID STAT QUEUE JOB_NAME SUBMIT_TIME",
             "-noheader"],
            capture_output=True, text=True, timeout=10)
        if cp.returncode != 0: return []
        jobs = []
        for line in cp.stdout.splitlines():
            parts = line.split(None, 4)
            if len(parts) >= 4:
                jobs.append({"jobid": parts[0], "stat": parts[1],
                             "queue": parts[2], "name": parts[3],
                             "submitted": parts[4] if len(parts)>4 else ""})
        return jobs
    except Exception:
        return []

def tail_file(path, n=50):
    p = Path(path)
    if not p.exists(): return ""
    try:
        with open(p) as f:
            lines = f.readlines()
        return "".join(lines[-n:])
    except Exception as e:
        return f"<could not read: {e}>"

def parse_md_progress(prod_out_path):
    p = Path(prod_out_path)
    if not p.exists(): return None
    try:
        with open(p) as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 10000))
            text = f.read()
        matches = re.findall(r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([\d.]+)", text)
        if matches:
            return int(matches[-1][0]), float(matches[-1][1])
    except Exception:
        pass
    return None

def sanitized_openfe_env(openfe_bin="openfe", base_env=None):
    """Return an environment dict in which `openfe` runs against its OWN conda
    env's toolchain, with any *system* Amber (e.g. /share/apps/.../amber*)
    stripped from PATH/PYTHONPATH and AMBERHOME pointed at the conda env.

    This fixes the "AmberTools version None" failure that happens when a
    detached launch inherits a login shell where `module load amber` /
    `source amber.sh` polluted the environment: the system antechamber shim
    then shadows the conda one even though the interactive shell looked fine.
    """
    import shutil
    env = dict(os.environ if base_env is None else base_env)

    ofe = shutil.which(openfe_bin, path=env.get("PATH")) or openfe_bin
    env_bin = Path(ofe).resolve().parent          # .../<env>/bin
    env_root = env_bin.parent                      # .../<env>

    def _is_system_amber_dir(p):
        if not p:
            return False
        try:
            rp = Path(p).resolve()
        except Exception:
            rp = Path(p)
        s = str(rp).lower()
        if str(env_root.resolve()).lower() in s:   # never strip our own env
            return False
        parts = set(rp.parts)
        amber_root = any(
            seg.startswith("amber") and "_" not in seg and "-" not in seg
            and (seg[5:].isdigit() or seg in ("amber", "ambertools"))
            for seg in parts)
        if not amber_root:
            return False
        return ("wrapped_progs" in parts) or ("bin" in parts) or \
               ("site-packages" in parts and "lib" in parts)

    def _strip_amber(value, sep=os.pathsep):
        return sep.join(p for p in (value or "").split(sep)
                        if p and not _is_system_amber_dir(p))

    # 1) ensure the conda env bin is FIRST on PATH, then drop system amber dirs
    path = env.get("PATH", "")
    path = _strip_amber(path)
    path = os.pathsep.join([str(env_bin)] +
                           [p for p in path.split(os.pathsep)
                            if p and p != str(env_bin)])
    env["PATH"] = path

    # 2) scrub PYTHONPATH of system amber (your trace imported parmed from
    #    /share/apps/.../amber22/lib/python3.11 -- a classic leak)
    if env.get("PYTHONPATH"):
        env["PYTHONPATH"] = _strip_amber(env["PYTHONPATH"])
        if not env["PYTHONPATH"]:
            env.pop("PYTHONPATH", None)

    # 3) point AMBERHOME at the conda env (antechamber lives there now)
    if (env_root / "bin" / "antechamber").exists():
        env["AMBERHOME"] = str(env_root)
    else:
        env.pop("AMBERHOME", None)

    return env


def spawn_detached(cmd, log_path, cwd=None, env=None):
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as f:
        f.write(f"\n# Spawned at {datetime.now().isoformat()}\n")
        f.write(f"# Command: {' '.join(str(c) for c in cmd)}\n\n")
        f.flush()
        proc = subprocess.Popen(cmd, cwd=cwd, env=env,
                                stdout=f, stderr=subprocess.STDOUT,
                                start_new_session=True)
    return proc.pid

def is_pid_alive(pid):
    if not pid: return False
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False

def workflow_root():
    here = Path(__file__).resolve()
    for parent in [here.parent.parent.parent, here.parent.parent]:
        if (parent / "run_amber.py").exists():
            return parent
    return Path.cwd()

def run_amber_py():
    return workflow_root() / "run_amber.py"

def prep_ligands_py():
    return workflow_root() / "prep_ligands.py"

def render_sidebar_globals():
    with st.sidebar:
        st.divider()
        st.subheader("Global defaults")
        proj = st.text_input("LSF project (-P)",
                             value=st.session_state.get("g_project","your-project"),
                             key="g_project")
        queue = st.text_input("GPU queue", value=st.session_state.get("g_queue","gpu"),
                              key="g_queue")
        walltime = st.text_input("Walltime (HH:MM)",
                                 value=st.session_state.get("g_walltime","24:00"),
                                 key="g_walltime")
        n_gpu = st.number_input("GPUs per job", 1, 4,
                                value=st.session_state.get("g_n_gpu",1), key="g_n_gpu")
    return {"project": proj, "queue": queue, "walltime": walltime, "n_gpu": int(n_gpu)}
