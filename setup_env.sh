#!/usr/bin/env bash
# Source this before running the pipeline:
#     source setup_env.sh
# Edit module names to match your cluster.

module purge

module load gcc/11.5  || echo "[setup_env] gcc/11.5 not available — adjust"
module load cuda      || echo "[setup_env] cuda module not available — adjust"
module load amber     || echo "[setup_env] amber module not available — adjust"

# Python >= 3.7 is required (stdlib 'dataclasses')
if ! python -c 'import sys; sys.exit(0 if sys.version_info >= (3,7) else 1)' 2>/dev/null; then
    module load python/3.11 2>/dev/null || \
    module load python/3.10 2>/dev/null || \
    module load python/3.9  2>/dev/null || \
    echo "[setup_env] No python>=3.7 module found — install miniconda."
fi

# Optional venv with extra Python packages
if [ -f "$HOME/venvs/amber-md/bin/activate" ]; then
    source "$HOME/venvs/amber-md/bin/activate"
fi

# Sanity
if [ -f "$(dirname "${BASH_SOURCE[0]}")/check_env.py" ]; then
    python "$(dirname "${BASH_SOURCE[0]}")/check_env.py" || true
fi
