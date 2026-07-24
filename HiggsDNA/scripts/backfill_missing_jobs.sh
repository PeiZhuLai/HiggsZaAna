#!/usr/bin/env bash
# Backfill missing chunks for HZa_merged Step B (data-ML friend parquet), then re-merge.
#
# "Missing" jobs come from RAL xrootd timeouts (retired) and watchdog-reaped stragglers.
# Per incomplete tag we re-run run_analysis.py: HiggsDNA's output-presence logic resubmits
# ONLY the jobs whose output_job_*_nominal.parquet is absent, reusing the ones already there.
# Reads are steered to FNAL (HDNA_REDIRECTOR) to dodge the RAL redirector that caused the
# original failures. After condor drains (with the same straggler-reap guard as the watchdog)
# we re-merge with merge_tag_outputs.py so merged_nominal reflects the new chunks.
#
# Usage:
#   backfill_missing_jobs.sh            # DRY-RUN: just print the reconcile backfill list
#   backfill_missing_jobs.sh --run [tag ...]
#                                       # actually backfill (all incomplete tags, or the given ones)
#
# Safe to run only after the main Step B (orchestrate_data_ml_production.sh) has finished
# all 32 tags -- do NOT run it concurrently with the driver (they would fight over condor).
set -uo pipefail

ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
PY=$ENVDIR/bin/python
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
export CONDA_PREFIX=$ENVDIR
export PATH=$ENVDIR/bin:/cvmfs/cms.cern.ch/common:/usr/bin:/bin
export PYTHONPATH=$REPO
export X509_USER_PROXY=${X509_USER_PROXY:-/tmp/x509up_u175325}
cd "$REPO"

FRIEND_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend
CFG_DIR=metadata/friend_tmp_configs
CATALOG=metadata/samples/za_merged_data_2024_perds.json
BF_LOG=${BF_LOG:-$REPO/logs/backfill_$(date +%Y%m%d).log}

REDIR=${REDIR:-root://cmsxrootd.fnal.gov}   # FNAL: avoids the flaky RAL redirector
FPO=${FPO:-2}
N_CORES=${N_CORES:-4}
POLL=${POLL:-45}
STABLE_SECONDS=${STABLE_SECONDS:-600}
STRAGGLER_MAX=${STRAGGLER_MAX:-15}

log(){ echo "[$(date +%m-%d\ %H:%M:%S)] $*" | tee -a "$BF_LOG"; }
# tag-scoped: count ONLY merged-data friend jobs (Iwd = .../HZa/HiggsZaAna/HiggsDNA),
# so co-running productions (HZgamma, Parquet2Rootfile, Sig_MC) never inflate the count
# or trip the "bulk done" straggler-reap. One merged tag runs at a time -> all such jobs
# belong to the current tag.
condor_total(){ condor_q pelai -af Iwd 2>/dev/null | grep -c "HZa/HiggsZaAna/HiggsDNA"; }

incomplete_tags(){
    $PY scripts/reconcile_stepB.py "$FRIEND_BASE" 2>/dev/null \
        | sed -n 's/^BACKFILL \([^ ]*\) .*/\1/p'
}

# ---- dry-run: show what would be backfilled ----
if [ "${1:-}" != "--run" ]; then
    echo "DRY-RUN. Incomplete tags (missing jobs) per reconcile_stepB.py:"
    $PY scripts/reconcile_stepB.py "$FRIEND_BASE" | sed -n '/tags=/,$p'
    echo
    echo "To actually backfill: $0 --run [tag ...]   (FNAL redirector: $REDIR)"
    exit 0
fi
shift

if [ -n "${1:-}" ]; then TAGS=("$@"); else mapfile -t TAGS < <(incomplete_tags); fi
if [ "${#TAGS[@]}" -eq 0 ]; then log "nothing to backfill -- all tags complete."; exit 0; fi

# refuse to run while the sequential driver is still going (condor contention)
if pgrep -f orchestrate_data_ml_production >/dev/null; then
    log "REFUSING: orchestrate_data_ml_production is still running. Let Step B finish first."
    exit 1
fi

log "backfill start: ${#TAGS[@]} tag(s) via REDIR=$REDIR -> ${TAGS[*]}"
export HDNA_REDIRECTOR="$REDIR"

for tag in "${TAGS[@]}"; do
    cfg="$CFG_DIR/${tag}.json"
    outdir="$FRIEND_BASE/${tag}"
    taskdir=$(ls -d "$outdir"/*_2024 2>/dev/null | head -1)
    if [ ! -f "$cfg" ] || [ -z "${taskdir:-}" ]; then
        log "[skip] $tag: missing config or task dir"; continue
    fi

    log "=== backfill $tag (FNAL) ==="
    # force DAS re-resolution with the FNAL redirector (drop cached infn URLs)
    for c in "${CATALOG%.json}"_sample_manager_full.json "$FRIEND_BASE/${tag}"/*_sample_manager_full.json; do
        [ -f "$c" ] && mv -f "$c" "$c.pre_backfill.$(date +%s)" 2>/dev/null && log "  cleared cache $c"
    done
    rm -f "$outdir/analysis_manager.pkl" "$outdir/analysis_manager_temp.pkl"

    # resubmit only the missing jobs (output-presence), monitor in background
    $PY scripts/run_analysis.py --config "$cfg" --log-level INFO --n_cores "$N_CORES" \
        --output_dir "$outdir" --batch_system condor --merge_outputs --fpo "$FPO" \
        >>"$BF_LOG" 2>&1 &
    ra_pid=$!

    # inline stall handling. CRITICAL: pre-existing chunks are ~days old, so we must NOT
    # judge "done" by their mtime. We (a) gate on `submitted` (condor was seen non-empty ->
    # run_analysis actually resubmitted the missing jobs), and (b) measure staleness only
    # against NEW production (chunks with mtime > t_start), else from t_start. A tag is done
    # when jobs were submitted, condor is drained to <=STRAGGLER_MAX, and no NEW chunk has
    # landed for STABLE_SECONDS (covers: new chunks settled, all-failed, or stuck stragglers).
    t_start=$(date +%s)
    max_wait=${MAX_WAIT_PER_TAG:-10800}   # 3h safety cap per tag
    submitted=0
    last_new=$t_start
    while kill -0 "$ra_pid" 2>/dev/null; do
        sleep "$POLL"
        now=$(date +%s)
        ct=$(condor_total); ct=${ct:-0}
        [ "$ct" != "0" ] && submitted=1
        newest=$(find "$taskdir" -name 'output_job_*_nominal.parquet' -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
        newest=${newest%.*}; newest=${newest:-0}
        if [ "$newest" -gt "$t_start" ] && [ "$newest" -gt "$last_new" ]; then last_new=$newest; submitted=1; fi

        # bail-out safety: run_analysis stuck without ever submitting
        if [ $((now - t_start)) -ge "$max_wait" ]; then
            log "  MAX_WAIT reached for $tag (submitted=$submitted) -> stop + remerge"
            kill -TERM "$ra_pid" 2>/dev/null; sleep 4; break
        fi
        [ "$submitted" = "0" ] && continue          # still resolving/submitting, wait
        [ "$ct" -gt "$STRAGGLER_MAX" ] && continue   # bulk still producing, wait

        idle=$((now - last_new))                     # secs since last NEW chunk (or since start)
        [ "$idle" -lt "$STABLE_SECONDS" ] && continue

        log "  bulk done ($ct job(s) left, ${idle}s since last new chunk) -> stop monitor + reap + remerge"
        kill -TERM "$ra_pid" 2>/dev/null
        sleep 4
        # reap ONLY this tag's leftover clusters (Iwd = .../HiggsDNA), never a blanket
        # 'condor_q pelai' which would also remove co-running HZgamma / Parquet2Rootfile jobs.
        leftover=$(condor_q pelai -af ClusterId Iwd 2>/dev/null | grep "HZa/HiggsZaAna/HiggsDNA" | awk '{print $1}' | sort -u | tr '\n' ' ')
        [ -n "${leftover// /}" ] && { log "  condor_rm $leftover"; condor_rm $leftover >>"$BF_LOG" 2>&1; }
        break
    done

    # always re-merge so merged_nominal reflects whatever chunks now exist
    if $PY scripts/merge_tag_outputs.py "$outdir" >>"$BF_LOG" 2>&1; then
        log "  re-merged $tag OK"
    else
        log "  re-merge FAILED for $tag"
    fi
done

log "backfill pass done. Reconcile:"
$PY scripts/reconcile_stepB.py "$FRIEND_BASE" | sed -n '/tags=/,$p' | tee -a "$BF_LOG"
log "Re-run '$0 --run' for another pass if tags remain incomplete (RAL may need retries)."
