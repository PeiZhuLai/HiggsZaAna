#!/usr/bin/env bash
# Start the CERN part of the HZa_merged lepclean re-production as systemd --user
# units, so it survives terminal-session restarts and ssh logouts (el9 lxplus has
# KillUserProcesses=yes; setsid nohup does not survive it).
#
# Units: hza-lc-w1..w4 (queue workers), hza-lc-sig (signal), hza-lc-watch (hold
# watchdog + heartbeat). Each starts with aklog: a unit is not a child of the ssh
# session and may otherwise have no AFS token (redirect to the AFS log fails and
# the unit exits 1 before writing anything).
#
# Dry-run by default; APPLY=1 to start.
#   bash launch.sh            # print what would start + pre-flight checks
#   APPLY=1 bash launch.sh
set -uo pipefail

HERE=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/prod_lepclean
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
N_WORKERS="${N_WORKERS:-4}"
APPLY="${APPLY:-0}"

# ---- pre-flight -----------------------------------------------------------
bad=0
tl=$(voms-proxy-info -file /tmp/x509up_u175325 -timeleft 2>/dev/null || echo 0)
echo "proxy /tmp: $((tl/3600)) h left"; [ "${tl}" -lt 172800 ] && { echo "  !! < 48 h"; bad=1; }
cmp -s /tmp/x509up_u175325 "${REPO}/x509up_u175325" || { echo "  !! repo proxy copy differs from /tmp"; bad=1; }
while read -r kind tag; do
    [ -z "${kind}" ] || [ "${kind:0:1}" = "#" ] && continue
    [ -f "${REPO}/metadata/parent_maps/${tag}.json" ] || { echo "  !! no parent map ${tag}"; bad=1; }
    [ -d "/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD/${tag}" ] || { echo "  !! no MLNanoAOD ${tag}"; bad=1; }
done < "${HERE}/tags_cern.txt"
for u in $(seq 1 "${N_WORKERS}") sig watch; do
    n="hza-lc-${u}"; [ "${u}" -eq "${u}" ] 2>/dev/null && n="hza-lc-w${u}"
    systemctl --user is-active --quiet "${n}" && { echo "  !! ${n} already active"; bad=1; }
done
echo "queue: $(grep -c '^data' "${HERE}/tags_cern.txt") data + $(grep -c '^bkg' "${HERE}/tags_cern.txt") bkg tags; claims so far: $(ls "${LOGD}/claims" 2>/dev/null | wc -l)"
[ "${bad}" -ne 0 ] && { echo "pre-flight FAILED"; exit 1; }
echo "pre-flight OK"

run() { echo "+ $*"; [ "${APPLY}" = "1" ] && "$@"; }
KRB="KRB5CCNAME=${KRB5CCNAME:-FILE:/run/user/$(id -u)/krb5cc}"
for i in $(seq 1 "${N_WORKERS}"); do
    run systemd-run --user --collect --unit="hza-lc-w${i}" --setenv="${KRB}" \
        /usr/bin/bash -c "aklog; exec /usr/bin/bash ${HERE}/queue_worker.sh ${i} >> ${LOGD}/worker_${i}.out 2>&1"
    [ "${APPLY}" = "1" ] && sleep 20   # stagger the first DAS queries / submissions
done
run systemd-run --user --collect --unit=hza-lc-sig --setenv="${KRB}" \
    /usr/bin/bash -c "aklog; exec /usr/bin/bash ${HERE}/run_signal.sh >> ${LOGD}/signal.out 2>&1"
run systemd-run --user --collect --unit=hza-lc-watch --setenv="${KRB}" --setenv=N_WORKERS="${N_WORKERS}" \
    /usr/bin/bash -c "aklog; exec /usr/bin/bash ${HERE}/hold_watchdog.sh >> ${LOGD}/watchdog.out 2>&1"
[ "${APPLY}" = "1" ] || echo "(dry-run; APPLY=1 to start)"
