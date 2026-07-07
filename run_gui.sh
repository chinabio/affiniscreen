#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

if [[ "${CONDA_DEFAULT_ENV:-}" != "amber-md" ]] && [[ -f "./activate_amber_md.sh" ]]; then
    source ./activate_amber_md.sh
fi

PORT="${1:-${PORT:-8501}}"
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"

echo "============================================================"
echo " AffiniScreen GUI starting on port $PORT"
echo " Open your browser to: http://localhost:$PORT"
echo "============================================================"

exec streamlit run "$(pwd)/amber_md/gui/Home.py" \
    --server.port "$PORT" \
    --server.address localhost \
    --browser.gatherUsageStats false \
    --server.headless true
