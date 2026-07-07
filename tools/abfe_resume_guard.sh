#!/usr/bin/env bash
# =============================================================================
# abfe_resume_guard.sh  (v3 / final47)  -- SUBMISSION-AWARE resume babysitter
#
# v1/v2 bug: each sweep ran `fep_driver --resume --submit`, which ALWAYS creates
# a new LSF array for incomplete windows. On a busy queue where nothing had
# finished yet, every 15-min sweep stacked another full ~81-job pipeline ->
# hundreds of duplicate jobs (observed: 889 = ~11 stacked submissions).
#
# v3 fix: BEFORE resubmitting, check whether THIS run already has live
# (PEND/RUN) fep_* jobs. If so, SKIP the sweep -- the windows are already
# queued, just not finished. Only resubmit when the queue has genuinely drained
# (preemption/maintenance killed jobs) but windows remain incomplete.
#
# It identifies "this run's" jobs two ways:
#   (a) the LSF array IDs recorded in <workdir>/fep/job_ids.json (written by the
#       driver each submit), and
#   (b) by job-name prefix fep_  owned by $USER  (fallback / belt-and-braces).
#
# USAGE
#   FIRST_LAUNCH=1 LIG_RESNAME=LIG nohup ./abfe_resume_guard.sh > guard.log 2>&1 &
#   ./abfe_resume_guard.sh /abs/path/to/workdir
# =============================================================================
set -uo pipefail

PROTEIN="${PROTEIN:-$HOME/Run_dir/protein.pdb}"
LIGAND="${LIGAND:-$HOME/Run_dir/ligand.mol2}"
LIG_RESNAME="${LIG_RESNAME:-LIG}"
LIG_CHARGE="${LIG_CHARGE:-0}"
CHARGE_METHOD="${CHARGE_METHOD:-bcc}"
PROJECT="${PROJECT:-your-project}"
QUEUE="${QUEUE:-gpu}"
WALLTIME="${WALLTIME:-48:00}"
NSTLIM_EQ="${NSTLIM_EQ:-250000}"
NSTLIM_PROD="${NSTLIM_PROD:-2500000}"
TEMP="${TEMP:-298.0}"
ACTIVATE="${ACTIVATE:-$HOME/Tools/amber_md_workflow_v2.5.0/activate_amber_md.sh}"

POLL_MIN="${POLL_MIN:-15}"
MAX_SWEEPS="${MAX_SWEEPS:-400}"
STALL_LIMIT="${STALL_LIMIT:-8}"
FIRST_SWEEP_DELAY_MIN="${FIRST_SWEEP_DELAY_MIN:-30}"

# shellcheck disable=SC1090
source "$ACTIVATE"

# ---- count THIS run's live (PEND|RUN) jobs --------------------------------
live_jobs() {
    # Prefer the array IDs the driver recorded for this workdir.
    local ids n=0 jid st
    local jfile="$WORKDIR/fep/job_ids.json"
    if [[ -s "$jfile" ]]; then
        ids=$(grep -oE '[0-9]{5,}' "$jfile" | sort -u)
        if [[ -n "$ids" ]]; then
            for jid in $ids; do
                # count PEND/RUN elements of this array id
                n=$((n + $(bjobs -w "$jid" 2>/dev/null | awk 'NR>1 && ($3=="PEND"||$3=="RUN")' | wc -l)))
            done
            echo "$n"; return
        fi
    fi
    # Fallback: any of MY fep_* jobs that are PEND/RUN.
    n=$(bjobs -w -J 'fep_*' 2>/dev/null | awk 'NR>1 && ($3=="PEND"||$3=="RUN")' | wc -l)
    echo "$n"
}

count_incomplete() {
    local n=0 po d leg
    for leg in complex_decharge complex_vdw solvent_decharge solvent_vdw; do
        for d in "$WORKDIR"/fep/"$leg"/lambda_*/; do
            [[ -d "$d" ]] || continue
            po="$d/prod.out"
            if [[ ! -s "$po" ]] || ! grep -qiE "Final Performance|wallclock|TIMINGS|Total wall time" "$po" 2>/dev/null; then
                n=$((n+1))
            fi
        done
    done
    echo "$n"
}

# ---- launch / resolve workdir ---------------------------------------------
if [[ "${FIRST_LAUNCH:-0}" == "1" ]]; then
    WORKDIR="$HOME/abfe_production_$(date +%Y%m%d_%H%M%S)"
    echo "[guard] FIRST_LAUNCH -> $WORKDIR (resname=$LIG_RESNAME)"
    python -m amber_md.fep_driver --mode abfe --work-dir "$WORKDIR" \
        --protein-pdb "$PROTEIN" --ligand-file "$LIGAND" \
        --ligand-resname "$LIG_RESNAME" --ligand-charge "$LIG_CHARGE" \
        --charge-method "$CHARGE_METHOD" --auto-boresch \
        --nstlim-eq "$NSTLIM_EQ" --nstlim-prod "$NSTLIM_PROD" --temperature "$TEMP" \
        --project "$PROJECT" --queue "$QUEUE" --walltime "$WALLTIME" --n-gpu 1 \
        --submit --analyze
    rc=$?
    [[ $rc -ne 0 ]] && { echo "[guard] launch failed rc=$rc"; exit $rc; }
    echo "[guard] launched; waiting ${FIRST_SWEEP_DELAY_MIN} min before first sweep ..."
    sleep "$((FIRST_SWEEP_DELAY_MIN*60))"
else
    WORKDIR="${1:?Usage: $0 /abs/path/to/workdir   (or FIRST_LAUNCH=1 $0)}"
fi

RESULT="$WORKDIR/fep/ABFE_RESULT.json"
echo "[guard] watching $WORKDIR (poll=${POLL_MIN}min stall=$STALL_LIMIT)"

prev=-1; stall=0; sweep=1
while [[ $sweep -le $MAX_SWEEPS ]]; do
    if [[ -s "$RESULT" ]]; then
        echo "[guard] ABFE_RESULT.json present -> QC and exit."
        python -m amber_md.abfe_qc "$WORKDIR" || true
        exit 0
    fi

    inc="$(count_incomplete)"
    alive="$(live_jobs)"
    echo "[guard] sweep $sweep : incomplete=$inc  live(PEND/RUN)=$alive"

    # *** THE FIX: never resubmit while this run still has live jobs. ***
    if [[ "$alive" -gt 0 ]]; then
        echo "[guard] $alive job(s) still queued/running -> NOT resubmitting; waiting."
        stall=0; prev="$inc"
        sleep "$((POLL_MIN*60))"; sweep=$((sweep+1)); continue
    fi

    # Queue drained but windows remain -> genuine preemption: stall accounting.
    if [[ "$inc" -eq "$prev" ]]; then stall=$((stall+1)); else stall=0; fi
    prev="$inc"
    if [[ "$stall" -ge "$STALL_LIMIT" ]]; then
        echo "[guard] queue empty but no progress for $stall sweeps -> real crash. Scan:"
        for leg in complex_decharge complex_vdw solvent_decharge solvent_vdw; do
            for d in "$WORKDIR"/fep/"$leg"/lambda_*/; do
                po="$d/prod.out"; [[ -s "$po" ]] || continue
                grep -qiE "NaN|vlimit|STOP PMEMD|\*\*\*\*|0 atoms" "$po" 2>/dev/null && {
                    echo "  !! $po"; grep -iE "NaN|vlimit|STOP PMEMD|\*\*\*\*|0 atoms" "$po" | tail -3; }
            done
        done
        exit 2
    fi

    echo "[guard] queue drained, $inc window(s) incomplete -> ONE resume submit."
    python -m amber_md.fep_driver --mode abfe --work-dir "$WORKDIR" \
        --ligand-file "$LIGAND" --ligand-resname "$LIG_RESNAME" \
        --project "$PROJECT" --queue "$QUEUE" --walltime "$WALLTIME" --n-gpu 1 \
        --resume --submit --analyze
    rc=$?
    [[ $rc -ne 0 ]] && echo "[guard] WARNING: --resume rc=$rc (continuing)"
    sleep "$((POLL_MIN*60))"; sweep=$((sweep+1))
done
echo "[guard] hit MAX_SWEEPS=$MAX_SWEEPS. Inspect $WORKDIR."
exit 3
