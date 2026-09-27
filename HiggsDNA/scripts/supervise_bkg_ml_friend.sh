#!/usr/bin/env bash
# Detached supervisor for the merged-ML friend production.
#
# Claude Code's background processes are children of its shell, so restarting or
# updating CC kills them. Condor jobs already in the queue survive (they run on
# worker nodes) but the driver that monitors and merges them does not -- jobs
# finish, outputs land, and nothing ever merges them.
#
# Launch with setsid so this outlives the session:
#   setsid nohup bash scripts/supervise_bkg_ml_friend.sh > /dev/null 2>&1 &
#
# What it does, every CHECK_INTERVAL:
#   * if a driver is running, leave it alone;
#   * if not, and tags are still unfinished, restart the driver for those tags
#     (run_merged_bkg_ml_friend.sh skips tags whose merged parquet exists, and
#     HiggsDNA skips jobs whose output already exists -- verified: the 2026-08-29
#     resume submitted 448 of 740 jobs, not all 740);
#   * when all four tags reconcile, write the DONE marker and stop.
#
# It does NOT run dataVmc -- plotting needs judgement, not automation.

set -uo pipefail

REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged
OUT="${BASE}/parquet_friend_ML"
STATE="${BASE}/SUPERVISOR_STATE.log"
DONE_MARK="${BASE}/SUPERVISOR_DONE"
DRIVER_LOG="${BASE}/bkg_ml_friend_supervised.log"
CHECK_INTERVAL="${CHECK_INTERVAL:-600}"
MAX_RESTARTS="${MAX_RESTARTS:-8}"
TAGS_ALL="Bkg_DYGto2LG_10to100_2024 Bkg_DYJetsTo2E_2024 Bkg_DYJetsTo2Mu_2024 Bkg_DYJetsTo2Tau_2024"

export PATH=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin:$PATH

log() { echo "$(date -Is) $*" >> "${STATE}"; }

# /tmp is wiped by reboots and by CC restarts, taking the grid proxy with it.
# Without it SampleManager.get_samples() dies on a bare `raise RuntimeError()`
# with no message -- the 2026-08-29 22:56 failure. A long-lived copy lives on
# AFS; restore from it rather than dying, and only give up if that is stale too.
PROXY_BACKUP="${PROXY_BACKUP:-/afs/cern.ch/user/p/pelai/x509up_condor}"
PROXY_PATH="/tmp/x509up_u$(id -u)"
ensure_proxy() {
    local left
    left=$(X509_USER_PROXY="${PROXY_PATH}" voms-proxy-info -timeleft 2>/dev/null || echo 0)
    if [ "${left:-0}" -gt 3600 ] 2>/dev/null; then
        return 0
    fi
    local bleft
    bleft=$(X509_USER_PROXY="${PROXY_BACKUP}" voms-proxy-info -timeleft 2>/dev/null || echo 0)
    if [ "${bleft:-0}" -gt 3600 ] 2>/dev/null; then
        cp -p "${PROXY_BACKUP}" "${PROXY_PATH}" && chmod 600 "${PROXY_PATH}"
        log "proxy restored from ${PROXY_BACKUP} (${bleft}s left)"
        return 0
    fi
    log "NO VALID PROXY: ${PROXY_PATH} has ${left:-0}s, backup has ${bleft:-0}s -- needs 'voms-proxy-init --rfc --voms cms -valid 192:00' from a human"
    return 1
}

merged_of() {  # $1=tag -> path or empty
    compgen -G "${OUT}/$1/*/merged_nominal.parquet" 2>/dev/null | head -1
}

pending_tags() {
    local out=""
    for t in ${TAGS_ALL}; do
        [ -z "$(merged_of "$t")" ] && out="${out}${t},"
    done
    echo "${out%,}"
}

driver_running() {
    # match the driver script, never this supervisor (pgrep -f would match both,
    # and matches the querying shell too -- check argv explicitly)
    pgrep -f "run_merged_bkg_ml_friend.sh" 2>/dev/null | grep -v "^$$\$" | head -1
}

log "=== supervisor start (pid $$, interval ${CHECK_INTERVAL}s) ==="
restarts=0

while true; do
    pend="$(pending_tags)"

    if [ -z "${pend}" ]; then
        log "all four tags have a merged parquet -- running reconciliation"
        rec="${BASE}/reconcile_final.txt"
        ( source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc12-opt/setup.sh 2>/dev/null
          python "${REPO}/scripts/reconcile_bkg_ml_friend.py" ) > "${rec}" 2>&1
        rc=$?
        log "reconciliation rc=${rc} -> ${rec}"
        tail -3 "${rec}" >> "${STATE}"
        if [ "${rc}" = "0" ]; then
            date -Is > "${DONE_MARK}"
            log "DONE -- all tags reconciled. Next: add_merged_flag.py then dataVmc."
            break
        fi
        log "reconciliation FAILED -- stopping for a human to look at ${rec}"
        break
    fi

    if ! ensure_proxy; then
        log "waiting for a valid proxy before touching the driver"
        sleep "${CHECK_INTERVAL}"
        continue
    fi

    drv="$(driver_running)"
    if [ -n "${drv}" ]; then
        now=$(date +%s)
        lm=$(stat -c %Y "${DRIVER_LOG}" 2>/dev/null || echo "${now}")
        idle=$(( (now - lm) / 60 ))
        log "driver ${drv} alive, log_idle=${idle}m, pending=${pend}"
        # HiggsDNA's run_condor_submit has no timeout on condor_submit; a wedged
        # schedd hangs it forever with every other signal looking healthy.
        if [ "${idle}" -ge 60 ]; then
            log "driver wedged (${idle}m idle) -- killing so it can be restarted"
            kill ${drv} 2>/dev/null
            pkill -P ${drv} 2>/dev/null
            sleep 20
        fi
    else
        if [ "${restarts}" -ge "${MAX_RESTARTS}" ]; then
            log "no driver and restart budget exhausted (${restarts}) -- stopping, needs a human"
            break
        fi
        restarts=$((restarts + 1))
        log "no driver running; restart ${restarts}/${MAX_RESTARTS} for tags: ${pend}"
        (
            cd "${REPO}" || exit 1
            TAGS="${pend}" FPO=4 N_CORES=4 BATCH_SYSTEM=condor \
                bash scripts/run_merged_bkg_ml_friend.sh
        ) >> "${DRIVER_LOG}" 2>&1 &
        sleep 60
    fi

    sleep "${CHECK_INTERVAL}"
done

log "=== supervisor exit ==="
