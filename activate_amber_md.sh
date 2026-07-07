#!/usr/bin/env bash
# Source this in every new shell:
#     source activate_amber_md.sh
#
# v2.5.2 changes:
#   - Prepend the workflow root (directory containing this script) to
#     PYTHONPATH so `import amber_md` works from any cwd, on the login
#     node or compute nodes. Previously the package was only importable
#     when you happened to be sitting in the workflow root, which broke
#     the in-LSF MM/GBSA report generator.
#
# v2.2.4 fixes:
#   - If we are already in 'amber-md', deactivate first so a fresh `module
#     load amber/22.8` sets AMBERHOME cleanly.
#   - Tolerate unset vars in the final echo.

# v2.5.2: resolve the workflow root from this script's own location.
# Works whether the script is sourced via absolute path, relative path,
# or a symlink. Falls back to PWD if BASH_SOURCE isn't populated (some
# minimal shells); that fallback matches the old, cwd-dependent behavior.
if [ -n "${BASH_SOURCE:-}" ]; then
    _AMBER_WORKFLOW_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
else
    _AMBER_WORKFLOW_ROOT="$PWD"
fi

# If we are already inside amber-md, drop back to base/none first so the
# module reload starts from a known-clean state.
if [ "${CONDA_DEFAULT_ENV:-}" = "amber-md" ]; then
    conda deactivate 2>/dev/null || true
fi

module purge
module load gcc/11.5
module load cuda/11.8           # adjust if pmemd.cuda complains
module load amber/22.8
module load pymol/3.0.4 2>/dev/null || true

# Save Amber-set environment before conda activate touches it
_AMBER_PYTHONPATH="${PYTHONPATH-}"
_AMBER_HOME="${AMBERHOME-}"
_AMBER_LDLIB="${LD_LIBRARY_PATH-}"
_AMBER_PERL5LIB="${PERL5LIB-}"

# Sanity warn if module load didn't set AMBERHOME
if [ -z "$_AMBER_HOME" ]; then
    echo "[activate] WARN: AMBERHOME not set after module load amber/22.8" >&2
    echo "[activate]       Try: module purge && module load amber/22.8" >&2
fi

# Initialize miniforge
if [ -f "$HOME/start_miniforge.sh" ]; then
    source "$HOME/start_miniforge.sh"
elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
else
    echo "[activate] ERROR: miniforge not found at $HOME/miniforge3" >&2
    return 1 2>/dev/null || exit 1
fi

# If base is already active, deactivate it first so we do not stack envs
if [ "${CONDA_DEFAULT_ENV:-}" = "base" ]; then
    conda deactivate
fi

conda activate amber-md

# Restore Amber-set environment that conda may have cleared
if [ -n "$_AMBER_PYTHONPATH" ]; then
    export PYTHONPATH="$_AMBER_PYTHONPATH${PYTHONPATH:+:$PYTHONPATH}"
fi
[ -n "$_AMBER_HOME" ]      && export AMBERHOME="$_AMBER_HOME"
[ -n "$_AMBER_PERL5LIB" ]  && export PERL5LIB="$_AMBER_PERL5LIB"
if [ -n "$_AMBER_LDLIB" ]; then
    export LD_LIBRARY_PATH="$_AMBER_LDLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

# v2.5.2: prepend the workflow root to PYTHONPATH so `import amber_md`
# works from any directory. Done AFTER restoring the Amber-set value, so
# the workflow root takes precedence (the amber_md package wins over
# anything in Amber's own python tree with the same name).
# Idempotent: if the root is already first on PYTHONPATH, do nothing.
if [ -d "$_AMBER_WORKFLOW_ROOT/amber_md" ]; then
    case ":${PYTHONPATH:-}:" in
        ":$_AMBER_WORKFLOW_ROOT:"*) ;;  # already first; no-op
        *) export PYTHONPATH="$_AMBER_WORKFLOW_ROOT${PYTHONPATH:+:$PYTHONPATH}" ;;
    esac
else
    echo "[activate] WARN: no amber_md/ subdir at $_AMBER_WORKFLOW_ROOT" >&2
    echo "[activate]       PYTHONPATH not extended; \`import amber_md\` may fail" >&2
fi

unset _AMBER_PYTHONPATH _AMBER_HOME _AMBER_LDLIB _AMBER_PERL5LIB _AMBER_WORKFLOW_ROOT

# Re-prepend AMBERHOME/bin to PATH so amber binaries beat conda's
if [ -n "${AMBERHOME:-}" ] && [ -d "$AMBERHOME/bin" ]; then
    case ":$PATH:" in
        *":$AMBERHOME/bin:"*) ;;
        *) export PATH="$AMBERHOME/bin:$PATH" ;;
    esac
fi

echo "[amber-md] env active"
echo "  python      = $(which python 2>/dev/null || echo MISSING)"
echo "  pmemd.cuda  = $(which pmemd.cuda 2>/dev/null || echo MISSING)"
echo "  tleap       = $(which tleap 2>/dev/null || echo MISSING)"
echo "  AMBERHOME   = ${AMBERHOME:-unset}"
echo "  conda env   = ${CONDA_DEFAULT_ENV:-none}"
echo "  PYTHONPATH  = ${PYTHONPATH:-unset}"
