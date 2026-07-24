#!/usr/bin/env bash
# Watchdog for HZa_merged Step B (orchestrate_data_ml_production.sh).
#
# The HiggsDNA condor monitor for a data task sometimes fails to recognise that
# all jobs are done and hangs in its poll loop -> never merges -> the sequential
# 32-tag driver stalls on that tag. This watchdog detects such a stall and:
#   1. faithfully merges that tag's completed chunks -> merged_nominal.parquet
#      (scripts/merge_tag_outputs.py, byte-equivalent to Task.merge_outputs for data)
#   2. kills that tag's run_analysis.py so the driver advances to the next tag.
#
# Stall condition (all must hold): the currently-running tag has NO merged_nominal
# yet, condor_q for pelai is empty (0 jobs total -> excludes held), and no new
# output chunk has been written for >= STABLE_SECONDS.
#
# Retired (RAL-timeout) jobs => missing chunks; we merge whatever completed and
# leave backfill (resume_missing_data.sh + FNAL redirector) for the reconciliation
# pass after all 32 tags are merged.
set -uo pipefail

ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
PY=$ENVDIR/bin/python
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
export CONDA_PREFIX=$ENVDIR
export PATH=$ENVDIR/bin:/usr/bin:/bin
export PYTHONPATH=$REPO
cd "$REPO"

FRIEND_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend
POLL=${POLL:-60}
STABLE_SECONDS=${STABLE_SECONDS:-600}   # no new chunk for this long => bulk production done
STRAGGLER_MAX=${STRAGGLER_MAX:-15}      # only reap when FEW jobs remain (true RAL stragglers);
                                        # many jobs left => tag still producing, keep waiting
STUCK_MINUTES=${STUCK_MINUTES:-75}      # (unused since reap-by-count) a running job older than this is RAL-stuck
WD_LOG=${WD_LOG:-$REPO/logs/watchdog_20260701.log}
log(){ echo "[$(date +%m-%d\ %H:%M:%S)] $*" | tee -a "$WD_LOG"; }

log "watchdog start (POLL=${POLL}s STABLE=${STABLE_SECONDS}s STUCK=${STUCK_MINUTES}min)"

# tag-scoped: count ONLY merged-data friend jobs (Iwd = .../HZa/HiggsZaAna/HiggsDNA), so
# co-running productions (HZgamma, Parquet2Rootfile, Sig_MC) don't block STALL detection.
condor_total(){ condor_q pelai -af Iwd 2>/dev/null | grep -c "HZa/HiggsZaAna/HiggsDNA"; }
# running clusters older than STUCK_MINUTES (RAL-stuck stragglers)
stuck_clusters(){
    condor_q pelai -run -af ClusterId JobStartDate ServerTime 2>/dev/null \
        | awk -v s="$STUCK_MINUTES" '($3-$2)/60>=s{print $1}'
}
merged_count(){ find "$FRIEND_BASE"/Data_2024_*/ -name merged_nominal.parquet 2>/dev/null | wc -l; }

while pgrep -f orchestrate_data_ml_production >/dev/null; do
    # current tag = output_dir of the live run_analysis.py
    # ONLY the merged-data run_analysis (parquet_friend/Data_2024_*). A concurrent unrelated
    # run_analysis (e.g. Sig_MC in parquet_DNA_tmp) must NOT be picked up or interfered with.
    outdir=$(ps -ef | grep "[r]un_analysis.py" | grep -oE -- "--output_dir [^ ]+" | awk '{print $2}' | grep "parquet_friend/Data_2024_" | head -1)
    if [ -z "${outdir:-}" ]; then
        sleep "$POLL"; continue
    fi
    tag=$(basename "$outdir")
    taskdir=$(ls -d "$outdir"/*_2024 2>/dev/null | head -1)
    [ -z "${taskdir:-}" ] && { sleep "$POLL"; continue; }

    # already merged? driver will advance on its own -> don't touch
    if [ -f "$taskdir/merged_nominal.parquet" ]; then
        sleep "$POLL"; continue
    fi

    # newest chunk age: if a chunk was written recently, production is still moving
    newest=$(find "$taskdir" -name 'output_job_*_nominal.parquet' -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
    if [ -z "${newest:-}" ]; then
        sleep "$POLL"; continue
    fi
    now=$(date +%s)
    age=$(awk -v n="$newest" -v now="$now" 'BEGIN{printf "%d", now-n}')

    # Production still moving? (new chunk within STABLE) -> wait.
    if [ "$age" -lt "$STABLE_SECONDS" ]; then
        sleep "$POLL"; continue
    fi

    # Bulk done: no new chunk for >=STABLE and no merged_nominal. Whatever remains is
    # not first-pass production (that writes chunks) but either a hung monitor (condor
    # empty) or RAL straggler retry churn (idle/running jobs that keep failing+resubmitting
    # and would each burn ~2h). Both advance the same way: merge what's done, move on.
    ct=$(condor_total); ct=${ct:-0}
    nstuck=$(stuck_clusters | grep -c .)          # running >= STUCK_MINUTES (RAL-hung)
    if [ "$ct" = "0" ]; then
        reason="condor-empty (monitor hung)"      # all jobs done, safe to merge
    elif [ "$ct" -le "$STRAGGLER_MAX" ]; then
        reason="$ct RAL straggler(s) left"        # few left => genuine stragglers, reap
    elif [ "$ct" -gt 0 ] && [ "$nstuck" -ge "$ct" ]; then
        # >STRAGGLER_MAX left BUT every one is a running job stuck >=STUCK_MINUTES (no young
        # idle jobs pending) => all RAL-hung, would each burn ~2h; reap. (Young idle jobs make
        # nstuck<ct, so a still-producing tag like the RunG false-reap case will NOT trigger.)
        reason="all $ct job(s) stuck >=${STUCK_MINUTES}min (RAL)"
    else
        # many jobs still queued/running (some young/idle): tag still producing, wait
        sleep "$POLL"; continue
    fi

    # ---- STALL: merge completed chunks, kill monitor (stops resubmits), reap leftovers ----
    log "STALL on $tag [$reason] chunk age ${age}s. Merging completed chunks."
    if $PY scripts/merge_tag_outputs.py "$outdir" >>"$WD_LOG" 2>&1; then
        pid=$(pgrep -f "run_analysis.py.*${tag}" | head -1)
        log "merge OK for $tag -> kill monitor pid ${pid:-none}; condor_rm leftover jobs; advance"
        [ -n "${pid:-}" ] && kill -TERM "$pid" 2>/dev/null
        sleep 4
        # monitor dead -> queued jobs won't resubmit; remove exactly the leftover cluster
        # ids of THIS tag (snapshot now), never a blanket 'condor_rm pelai'.
        # remove ONLY this tag's leftover jobs (match Iwd to the tag), never a blanket
        # pelai sweep -- a concurrent Sig_MC production shares the queue and must be spared.
        leftover=$(condor_q pelai -af ClusterId Iwd 2>/dev/null | grep "$tag" | awk '{print $1}' | sort -u | tr '\n' ' ')
        if [ -n "${leftover// /}" ]; then
            log "condor_rm leftover clusters: $leftover"
            condor_rm $leftover >>"$WD_LOG" 2>&1
        fi
        sleep 6
        log "merged_nominal now $(merged_count)/32"
    else
        log "MERGE FAILED for $tag -- leaving monitor alive, will retry next cycle"
    fi
    sleep "$POLL"
done

log "orchestrator exited. final merged_nominal $(merged_count)/32. watchdog stop."
