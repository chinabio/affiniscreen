#!/bin/bash
# ============================================================================
# submit_mmgbsa_for_screen.sh
# Submits one MM/GBSA LSF job per lig_*/ in a batch screen directory,
# chained to start AFTER each ligand's MD job finishes.
#
# Usage:
#   cd ~/Run_dir/Amber_run/screen_20260525_0932
#   bash submit_mmgbsa_for_screen.sh
#
# Env knobs:
#   QUEUE    (default gpu  ; set normalQ to offload to CPU nodes)
#   WALLTIME (default 2:00)
#   PROJECT  (default your-project)
#   NCPU     (default 1    ; set >1 to use MMPBSA.py.MPI)
#   SALTCON  (default 0.15)
#   IGB      (default 8)
#   ACTIVATE (default $HOME/Tools/affiniscreen/activate_amber_md.sh)
# ============================================================================
set -euo pipefail

QUEUE="${QUEUE:-gpu}"
WALLTIME="${WALLTIME:-2:00}"
PROJECT="${PROJECT:-your-project}"
NCPU="${NCPU:-1}"
SALTCON="${SALTCON:-0.15}"
IGB="${IGB:-8}"
ACTIVATE="${ACTIVATE:-$HOME/Tools/affiniscreen/activate_amber_md.sh}"

if [ "$NCPU" -gt 1 ]; then
    MMPBSA_CMD="mpirun -np ${NCPU} MMPBSA.py.MPI"
else
    MMPBSA_CMD="MMPBSA.py"
fi

shopt -s nullglob
LIGDIRS=( lig_*/ )
if [ ${#LIGDIRS[@]} -eq 0 ]; then
    echo "ERROR: no lig_*/ directories found in $(pwd)" >&2
    exit 1
fi

for d in "${LIGDIRS[@]}"; do
    wd=$(realpath "$d")
    name=$(basename "$d")
    echo "==> $name"

    if [ -s "$wd/mmgbsa/FINAL_RESULTS_MMPBSA.dat" ]; then
        echo "    already complete"
        continue
    fi

    for f in build/complex.prmtop topo/complex.prmtop topo/receptor.prmtop topo/ligand.prmtop; do
        if [ ! -s "$wd/$f" ]; then
            echo "    SKIP: missing $wd/$f" >&2
            continue 2
        fi
    done

    mkdir -p "$wd/mmgbsa"

    if [ ! -s "$wd/mmgbsa/mmgbsa.in" ]; then
        cat > "$wd/mmgbsa/mmgbsa.in" <<EOF
MM/GBSA (workaround script)
&general
  startframe=1, endframe=9999999, interval=1, verbose=2, keep_files=0
/
&gb
  igb=${IGB}, saltcon=${SALTCON}
/
EOF
    fi

    cat > "$wd/jobs/run_mmgbsa.lsf" <<EOF
#!/bin/bash
#BSUB -q ${QUEUE}
#BSUB -P ${PROJECT}
#BSUB -J MMGBSA_${name}
#BSUB -W ${WALLTIME}
#BSUB -n ${NCPU}
#BSUB -o mmgbsa.%J.out
#BSUB -e mmgbsa.%J.err

module purge
module load gcc/11.5
module load amber/22.8
source ${ACTIVATE}
set -euo pipefail

if [ ! -s ${wd}/jobs/prod.nc ]; then
    echo "[MM/GBSA] prod.nc missing - aborting." >&2
    exit 1
fi

mkdir -p ${wd}/mmgbsa
cd ${wd}/mmgbsa

cpptraj -p ${wd}/build/complex.prmtop <<'E'
trajin ${wd}/jobs/prod.nc
strip :WAT,Na+,Cl-,K+,Mg+2,Ca+2
autoimage
trajout prod_dry.nc netcdf
go
E

${MMPBSA_CMD} -O \\
  -i  ${wd}/mmgbsa/mmgbsa.in \\
  -o  ${wd}/mmgbsa/FINAL_RESULTS_MMPBSA.dat \\
  -sp ${wd}/build/complex.prmtop \\
  -cp ${wd}/topo/complex.prmtop \\
  -rp ${wd}/topo/receptor.prmtop \\
  -lp ${wd}/topo/ligand.prmtop \\
  -y  ${wd}/mmgbsa/prod_dry.nc

echo "[MM/GBSA] done \$(date)"
EOF

    chmod +x "$wd/jobs/run_mmgbsa.lsf"

    md_jid=""
    if [ -s "$wd/pipeline.log" ]; then
        md_jid=$(grep -oE 'Job <[0-9]+>' "$wd/pipeline.log" | head -1 | grep -oE '[0-9]+' || true)
    fi

    if [ -n "$md_jid" ]; then
        echo "    bsub -w done($md_jid) < jobs/run_mmgbsa.lsf"
        ( cd "$wd/jobs" && bsub -w "done($md_jid)" < run_mmgbsa.lsf )
    else
        echo "    WARN: no MD jobid found - submitting without dependency"
        ( cd "$wd/jobs" && bsub < run_mmgbsa.lsf )
    fi
done

echo
echo "All MM/GBSA jobs submitted. Monitor with:  bjobs -J 'MMGBSA_*'"
