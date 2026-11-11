#!/usr/bin/env bash
# Memory-hold handler for run_signal_ml.sh (one condor job per nano file, 3000 MB).
# Only jobs whose Cmd is under the HiggsDNA eos_logs of the signal samples are touched.
# Memory hold: RequestMemory 3000 -> 9000 -> 12000 for that job alone (qedit + release),
# and its batch_submit txt is raised too (only up) so a HiggsDNA resubmission keeps it.
# Any other hold, or a memory hold already at 12000: exit 3 so the caller looks at it.
# Exits 0 when the given systemd unit is no longer active, 4 at the time cap.
#   bash watch_sig_ml.sh <unit> [cap_seconds]
set -uo pipefail
UNIT="${1:?unit}"
CAP="${2:-1800}"
LOGS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/Sig_MC_MLNANO_all
ACT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean/watchdog_actions.log
end=$(( $(date +%s) + CAP ))
while systemctl --user is-active --quiet "${UNIT}"; do
    while IFS=$'\t' read -r id mem cmd reason; do
        [[ "${cmd}" == "${LOGS}/"* ]] || continue
        if [[ "${reason}" == *memory* ]] && [ "${mem}" -lt 12000 ]; then
            new=$([ "${mem}" -lt 9000 ] && echo 9000 || echo 12000)
            sub=$(ls "$(dirname "${cmd}")"/*_batch_submit_job*.txt 2>/dev/null | head -1)
            [ -n "${sub}" ] && sed -i -E "s/^(RequestMemory = )[0-9]+$/\1${new}/" "${sub}"
            condor_qedit "${id}" RequestMemory "${new}" >/dev/null && condor_release "${id}" >/dev/null
            echo "$(date '+%F %T') ${id} $(basename "$(dirname "$(dirname "${cmd}")")")/$(basename "$(dirname "${cmd}")") HOLD memory ${mem} -> ${new} + release :: ${reason:0:120}" | tee -a "${ACT}"
        else
            echo "UNHANDLED HOLD ${id} mem=${mem} :: ${reason}"; exit 3
        fi
    done < <(condor_q -cons 'JobStatus==5' -af:t ClusterId RequestMemory Cmd HoldReason 2>/dev/null)
    [ "$(date +%s)" -ge "${end}" ] && exit 4
    sleep 60
done
exit 0
