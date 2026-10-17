#!/usr/bin/env bash
# Hold watchdog + heartbeat for the HZa_merged lepclean re-production (CERN).
#
# HiggsDNA's submit template has OnExitHold = ExitStatus != 0, and CERN never
# releases held jobs by itself, while run_analysis waits on held jobs forever.
# Every POLL seconds, for held jobs of THIS production only (Cmd under the HZa
# HiggsDNA eos_logs AND batch name in this production's set):
#   * memory (HoldReason mentions memory, or MemoryUsage >= 90% of the request):
#       3000 -> 6000 -> 12000 -> 16000 MB via condor_qedit + condor_release, and
#       the per-job submit file is raised too (only up, never down) so a later
#       HiggsDNA resubmission keeps it. Released within one poll, so the held
#       time does not eat the walltime. Above 16000: left held, ALERT.
#   * anything else (exit != 0: transient xrootd/EOS, strict MLPhoton join):
#       release up to 3 times per job, then condor_rm so HiggsDNA resubmits it as
#       a new attempt (it retires a job after 5; the worker then unretires once).
# Each action is logged with job id, reason, old/new memory.
#
# Heartbeat line every poll in heartbeat.log; exits when every worker and the
# signal unit wrote their done marker and nothing of ours is left in the queue.
set -uo pipefail

HERE=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/prod_lepclean
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
ACT="${LOGD}/watchdog_actions.log"
HB="${LOGD}/heartbeat.log"
STATE="${LOGD}/release_counts.txt"
POLL="${POLL:-600}"
N_WORKERS="${N_WORKERS:-4}"
export PATH="/usr/bin:/bin:/usr/sbin:/sbin"
touch "${STATE}"

# batch names HiggsDNA gives this production's jobs: <sample>_<year>
NAMES="$(awk '$1=="data"{print $2"_2024"} $1=="bkg"{s=$2; sub(/^Bkg_/,"",s); print s}' "${HERE}/tags_cern.txt"
         for m in 1 2 3 4 5 6 7 8 9; do echo "mA_MLNANO_M0p${m}_2024"; done)"

while true; do
    q=$(timeout 120 condor_q -af:t ClusterId ProcId JobStatus JobBatchName RequestMemory MemoryUsage Cmd RemoteUserCpu JobCurrentStartDate ResidentSetSize_RAW HoldReason 2>/dev/null)
    if [ $? -ne 0 ]; then echo "$(date '+%F %T') condor_q failed" >> "${HB}"; sleep 60; continue; fi
    ours=$(echo "${q}" | awk -F'\t' -v names="${NAMES}" '
        BEGIN{n=split(names,a,"\n"); for(i=1;i<=n;i++) ok[a[i]]=1}
        $7 ~ /HiggsZaAna\/HiggsDNA\/eos_logs\// && ($4 in ok)')
    n_idle=$(echo "${ours}" | awk -F'\t' '$3==1' | grep -c . || true)
    n_run=$(echo "${ours}"  | awk -F'\t' '$3==2' | grep -c . || true)
    n_held=$(echo "${ours}" | awk -F'\t' '$3==5' | grep -c . || true)

    echo "${ours}" | awk -F'\t' '$3==5' | while IFS=$'\t' read -r cid pid st bn req use cmd cpu t0 rss reason; do
        id="${cid}.${pid}"
        use=${use%.*}; [ "${use}" = "undefined" ] && use=0
        sub="${cmd%_executable_job*}_batch_submit_job${cmd##*_executable_job}"; sub="${sub%.sh}.txt"
        if echo "${reason}" | grep -qi 'memory' || [ "${use:-0}" -ge $(( req * 9 / 10 )) ]; then
            # next rung above the request that also covers 1.3x the measured peak
            # (a job seen at 7325 MB must not be "raised" to 6000)
            want=$(( use * 13 / 10 )); new=0
            for rung in 6000 12000 16000; do
                if [ "${rung}" -gt "${req}" ] && [ "${rung}" -ge "${want}" ]; then new=${rung}; break; fi
            done
            [ "${new}" -eq 0 ] && [ "${req}" -lt 16000 ] && new=16000
            if [ "${new}" -eq 0 ]; then
                if ! grep -q "^ALERT ${id} " "${STATE}"; then
                    echo "ALERT ${id} " >> "${STATE}"
                    echo "$(date '+%F %T') ALERT ${id} ${bn} still memory-held at ${req} MB (use=${use}) -- left held" >> "${ACT}"
                    echo "$(date '+%F %H:%M') ALERT prod_lepclean ${id} ${bn} memory-held at ${req} MB" >> "${REG}"
                fi
                continue
            fi
            condor_qedit "${id}" RequestMemory "${new}" >/dev/null 2>&1
            if [ -f "${sub}" ]; then
                cur=$(awk '/^RequestMemory/{print $3}' "${sub}")
                [ "${cur:-0}" -lt "${new}" ] && sed -i "s/^RequestMemory = .*/RequestMemory = ${new}/" "${sub}"
            fi
            condor_release "${id}" >/dev/null 2>&1
            echo "$(date '+%F %T') MEM ${id} ${bn} req=${req} use=${use} -> ${new} | ${reason:0:120}" >> "${ACT}"
        else
            n=$(grep -c "^REL ${id}$" "${STATE}" || true)
            if [ "${n}" -lt 3 ]; then
                echo "REL ${id}" >> "${STATE}"
                condor_release "${id}" >/dev/null 2>&1
                echo "$(date '+%F %T') RELEASE#$((n+1)) ${id} ${bn} | ${reason:0:160}" >> "${ACT}"
            else
                condor_rm "${id}" >/dev/null 2>&1
                echo "$(date '+%F %T') RM ${id} ${bn} after 3 releases (HiggsDNA resubmits) | ${reason:0:160}" >> "${ACT}"
            fi
        fi
    done

    # stalled: running > STALL_MIN wall minutes with < 5 CPU minutes. Seen 2026-10-14:
    # xrdcp from xrootd-cms.infn.it hung 110 min (no timeout) and would sit until the
    # 8 h walltime. condor_rm -> HiggsDNA resubmits it as a new attempt.
    echo "${q}" | awk -F'\t' -v names="${NAMES}" -v now="$(date +%s)" -v lim="${STALL_MIN:-60}" '
        BEGIN{n=split(names,a,"\n"); for(i=1;i<=n;i++) ok[a[i]]=1}
        $7 ~ /HiggsZaAna\/HiggsDNA\/eos_logs\// && ($4 in ok) && $3==2 && $9 != "undefined" &&
        (now-$9)/60 >= lim && $8+0 < 300 {printf "%s.%s\t%s\t%d\t%d\n", $1, $2, $4, (now-$9)/60, $8}' |
    while IFS=$'\t' read -r id bn wall cpu; do
        condor_rm "${id}" >/dev/null 2>&1
        echo "$(date '+%F %T') STALL_RM ${id} ${bn} wall=${wall}min cpu=${cpu}s (HiggsDNA resubmits)" >> "${ACT}"
    done

    # running but over its memory request: the cgroup does not kill it, the job just
    # stops accumulating CPU (seen 2026-10-14: RSS 6 GB on a 3000 MB request, CPU frozen
    # for >2 h, never held). Same ladder as a memory hold, applied via hold/qedit/release.
    # Only when it is also stuck: >= 60 min wall and CPU efficiency < 0.3 (a job just
    # above its request is often still progressing and must not be restarted).
    echo "${q}" | awk -F'\t' -v names="${NAMES}" -v now="$(date +%s)" '
        BEGIN{n=split(names,a,"\n"); for(i=1;i<=n;i++) ok[a[i]]=1}
        $7 ~ /HiggsZaAna\/HiggsDNA\/eos_logs\// && ($4 in ok) && $3==2 && $10 != "undefined" &&
        $9 != "undefined" && $10/1024 > $5 && (now-$9) >= 3600 && $8 < 0.3*(now-$9) {printf "%s.%s\t%s\t%d\t%d\t%s\n", $1, $2, $4, $5, $10/1024, $7}' |
    while IFS=$'\t' read -r id bn req rssmb cmd; do
        want=$(( rssmb * 13 / 10 )); new=0
        for rung in 6000 12000 16000; do
            if [ "${rung}" -gt "${req}" ] && [ "${rung}" -ge "${want}" ]; then new=${rung}; break; fi
        done
        [ "${new}" -eq 0 ] && [ "${req}" -lt 16000 ] && new=16000
        [ "${new}" -eq 0 ] && continue
        sub="${cmd%_executable_job*}_batch_submit_job${cmd##*_executable_job}"; sub="${sub%.sh}.txt"
        condor_hold "${id}" >/dev/null 2>&1
        condor_qedit "${id}" RequestMemory "${new}" >/dev/null 2>&1
        if [ -f "${sub}" ]; then
            cur=$(awk '/^RequestMemory/{print $3}' "${sub}")
            [ "${cur:-0}" -lt "${new}" ] && sed -i "s/^RequestMemory = .*/RequestMemory = ${new}/" "${sub}"
        fi
        condor_release "${id}" >/dev/null 2>&1
        echo "$(date '+%F %T') OVERMEM ${id} ${bn} running rss=${rssmb}MB > req=${req} -> ${new} (hold/qedit/release)" >> "${ACT}"
    done

    nd=$(grep -c 'WORKER_DONE' "${LOGD}/status.txt" 2>/dev/null || true)
    sd=$(grep -c 'SIGNAL_DONE' "${LOGD}/status.txt" 2>/dev/null || true)
    np=$(grep -cE ' (DONE|FAIL) ' "${LOGD}/status.txt" 2>/dev/null || true)
    tl=$(voms-proxy-info -file /tmp/x509up_u175325 -timeleft 2>/dev/null || echo 0)
    echo "$(date '+%F %T') idle=${n_idle} run=${n_run} held=${n_held} tags_finished=${np} workers_done=${nd}/${N_WORKERS} signal_done=${sd} proxy_h=$((tl/3600))" >> "${HB}"
    if [ "${tl}" -lt 86400 ] && ! grep -q PROXY_LOW "${STATE}"; then
        echo "PROXY_LOW" >> "${STATE}"
        echo "$(date '+%F %H:%M') ALERT prod_lepclean grid proxy < 24h left -- voms-proxy-init needed" >> "${REG}"
    fi
    if [ "${nd}" -ge "${N_WORKERS}" ] && [ "${sd}" -ge 1 ] && [ $((n_idle + n_run + n_held)) -eq 0 ]; then
        echo "$(date '+%F %H:%M') DONE  prod_lepclean CERN part -> ${LOGD}/status.txt" >> "${REG}"
        echo "$(date '+%F %T') all done, exiting" >> "${HB}"
        exit 0
    fi
    sleep "${POLL}"
done
