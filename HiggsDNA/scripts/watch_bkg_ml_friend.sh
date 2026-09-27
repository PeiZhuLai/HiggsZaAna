#!/usr/bin/env bash
# Watchdog for the background merged-ML friend production.
#
# Records one status line per interval so the run stays auditable even if the
# interactive session dies. It does NOT act on its own -- run_analysis.py already
# re-queues unfinished jobs -- it records, and flags held jobs for a human.
#
# condor_q failing is NOT an empty queue: a loaded schedd returns rc!=0, and
# `condor_q | wc -l` would read that as 0 and declare the run finished. Every
# reading here checks rc first and records "QUERY-FAILED" instead of a count.
#
# Usage: nohup bash scripts/watch_bkg_ml_friend.sh <driver_pid> [interval_s] &

set -uo pipefail

DRIVER_PID="${1:?usage: watch_bkg_ml_friend.sh <driver_pid> [interval_s]}"
INTERVAL="${2:-600}"
OUT_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML
STATUS=/eos/cms/store/group/phys_susy/pelai/HZa_merged/bkg_ml_friend_watch.log
DRIVER_LOG="${DRIVER_LOG:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/bkg_ml_friend_full.log}"
STALE_MIN="${STALE_MIN:-30}"
TAGS="Bkg_DYGto2LG_10to100_2024 Bkg_DYJetsTo2E_2024 Bkg_DYJetsTo2Mu_2024 Bkg_DYJetsTo2Tau_2024"

echo "=== watchdog start $(date -Is) driver=${DRIVER_PID} interval=${INTERVAL}s ===" >> "${STATUS}"

while true; do
    ts=$(date -Is)

    if kill -0 "${DRIVER_PID}" 2>/dev/null; then
        drv="alive"
    else
        drv="GONE"
    fi

    q=$(condor_q -totals 2>/dev/null | grep "for pelai")
    if [ $? -ne 0 ] || [ -z "${q}" ]; then
        counts="QUERY-FAILED"
    else
        counts=$(echo "${q}" | sed 's/Total for pelai: //')
    fi

    held=$(condor_q -constraint 'JobStatus==5' -af ClusterId HoldReasonCode 2>/dev/null | wc -l)

    merged=""
    for t in ${TAGS}; do
        if [ -f "${OUT_BASE}/${t}/${t%_2024}_2024/merged_nominal.parquet" ] \
        || compgen -G "${OUT_BASE}/${t}/*/merged_nominal.parquet" > /dev/null 2>&1; then
            merged="${merged} ${t%%_2024}:DONE"
        else
            n=$(compgen -G "${OUT_BASE}/${t}/*/*.parquet" 2>/dev/null | wc -l)
            merged="${merged} ${t%%_2024}:${n}"
        fi
    done

    # Log staleness. A driver can sit in do_wait forever: HiggsDNA's
    # run_condor_submit calls do_cmd("condor_submit ...") with no timeout, so a
    # wedged schedd hangs it indefinitely. The process stays alive, the queue
    # stays empty, and every other signal here reads normal -- on 2026-08-29 that
    # went unnoticed for 9.5 hours. Only the log going quiet reveals it.
    now=$(date +%s)
    log_mtime=$(stat -c %Y "${DRIVER_LOG}" 2>/dev/null || echo "${now}")
    stale=$(( (now - log_mtime) / 60 ))

    echo "${ts} driver=${drv} held=${held} log_idle=${stale}m | ${counts} |${merged}" >> "${STATUS}"

    if [ "${stale}" -ge "${STALE_MIN}" ] 2>/dev/null; then
        echo "${ts}   STALE: driver log untouched for ${stale} min (threshold ${STALE_MIN}) -- likely wedged in condor_submit/condor_q" >> "${STATUS}"
    fi

    if [ "${held}" -gt 0 ] 2>/dev/null; then
        echo "${ts}   HELD DETAIL:" >> "${STATUS}"
        condor_q -constraint 'JobStatus==5' -af ClusterId ProcId HoldReason 2>/dev/null \
            | sort | uniq -c | sort -rn | head -5 | sed 's/^/    /' >> "${STATUS}"
    fi

    if [ "${drv}" = "GONE" ]; then
        echo "${ts} driver exited -- watchdog stopping" >> "${STATUS}"
        break
    fi

    sleep "${INTERVAL}"
done
