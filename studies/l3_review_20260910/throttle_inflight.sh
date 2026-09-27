#!/usr/bin/env bash
#
# Bound the number of HZa production jobs in flight, by holding the excess and
# releasing them as slots free up.
#
# WHY
#   Measured on 2026-09-13 at ~6100 concurrent jobs: THREE randomly sampled
#   running jobs all had rchar frozen over a 15 s window and sat in
#   futex_do_wait; two had produced 81 bytes of stdout after 58 minutes, i.e.
#   they were stuck inside the FIRST uproot.open. 6118 running jobs produced one
#   chunk in 30 minutes.
#   The same file opened from the login node, single client, in 8.6 s
#   (213122 events). So neither the file, nor the redirector, nor the proxy, nor
#   B-mode is at fault -- the failure is purely a function of how many clients
#   hit root://xrootd-cms.infn.it at once. fsspec's async loop has no working
#   timeout there (uproot's timeout=300 does not reach it), so a stalled open
#   never returns and the job burns its whole walltime doing nothing.
#   This also re-explains the 4515 jobs "removed for wall time" on 2026-09-12:
#   they were deadlocked the same way, not victims of the expired proxy. (The
#   2059 exit-1 [3010] jobs WERE the proxy -- their stack traces are stamped
#   after it expired.)
#
# HOW
#   Jobs are already submitted and configured, so nothing is resubmitted: the
#   excess is parked with condor_hold and let back in a few hundred at a time.
#   User holds get HoldReasonCode 1, which release_memory_holds.sh ignores (it
#   only touches code 34), so the two loops do not fight over the same jobs.
#
# USAGE
#   bash throttle_inflight.sh                        # dry run, report only
#   DRY_RUN=0 HOLD_ALL_FIRST=1 bash throttle_inflight.sh   # park everything, then meter
#   DRY_RUN=0 LOOP=1 nohup bash throttle_inflight.sh >/dev/null 2>&1 &
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
LOOP="${LOOP:-0}"
TARGET="${TARGET:-500}"
INTERVAL="${INTERVAL:-180}"
HOLD_ALL_FIRST="${HOLD_ALL_FIRST:-0}"
OUT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/throttle.log"

# Anchored, and Owner-scoped. The schedd is shared: an unanchored batch-name
# regex has already matched 1516 of another user's jobs once.
# Select by EXECUTABLE PATH, not by batch name. At 20:07 on 2026-09-13 an
# unrelated HZgamma production (DYto2L-2Jets_MLL-50_*J_FxFx_*) put 5400 jobs on
# the same schedd; the old regexp("^DY",JobBatchName) counted every one of them
# as "in flight", so room went permanently negative and this throttle stopped
# releasing its own parked jobs -- the HZa production sat at zero progress while
# looking busy. A HOLD_ALL_FIRST pass under that regex would have parked the
# other project's jobs outright.
# Cmd carries the repo path and cannot collide:
#   HZa : /afs/.../HZa/HiggsZaAna/HiggsDNA/eos_logs/...
#   HZg : /afs/.../HZgamma/higgsdna-hzg-run3/...
CONST='Owner=="pelai" && regexp("HZa/HiggsZaAna",Cmd)'

n_status() { condor_q -constraint "${CONST} && JobStatus==$1" -af ClusterId 2>/dev/null | wc -l | tr -d ' '; }

log() { echo "$(date '+%F %T') $*" >> "$OUT"; }

if [ "${HOLD_ALL_FIRST}" = "1" ]; then
    idle=$(n_status 1); run=$(n_status 2)
    log "parking everything: ${idle} idle + ${run} running (all deadlocked on xrootd)"
    if [ "$DRY_RUN" = "1" ]; then
        echo "would run: condor_hold -constraint '${CONST} && (JobStatus==1 || JobStatus==2)'"
    else
        condor_hold -constraint "${CONST} && (JobStatus==1 || JobStatus==2)" >> "$OUT" 2>&1
    fi
fi

meter() {
    local idle run held inflight room
    idle=$(n_status 1); run=$(n_status 2)
    # Only OUR parked jobs are releasable here. Code 34 (memory) belongs to
    # release_memory_holds.sh, which bumps RequestMemory before releasing --
    # releasing one here without the bump would just re-hold it immediately.
    held=$(condor_q -constraint "${CONST} && JobStatus==5 && HoldReasonCode==1" -af ClusterId ProcId 2>/dev/null \
           | awk '{print $1"."$2}')
    inflight=$((idle + run))
    local nheld; nheld=$(wc -w <<< "${held:-}")
    room=$((TARGET - inflight))
    log "inflight=${inflight} (idle=${idle} run=${run})  parked=${nheld}  room=${room}"
    [ -z "${held// }" ] && { log "nothing parked left to release"; return 1; }
    if [ "$room" -le 0 ]; then return 0; fi
    local batch; batch=$(tr ' ' '\n' <<< "$held" | grep -v '^$' | head -n "$room")
    local cnt; cnt=$(wc -w <<< "$batch")
    log "releasing ${cnt}"
    if [ "$DRY_RUN" = "1" ]; then
        echo "would run: condor_release <${cnt} ids>"
    else
        xargs -n 50 <<< "$batch" | while read -r b; do condor_release $b >> "$OUT" 2>&1; done
    fi
    return 0
}

while :; do
    meter || break
    [ "$LOOP" = "1" ] && [ "$DRY_RUN" != "1" ] || break
    # Deliberately NOT tied to the driver being alive. The drivers were killed on
    # purpose: HiggsDNA's retry logic calls condor_release on held jobs, so with a
    # driver running the parked pool drained itself (5577 -> 3104 parked while
    # in-flight climbed 500 -> 2952, which is how the first throttle attempt lost
    # control). Jobs already have their config and submit file and run fine with no
    # driver; merging is done by one final driver pass afterwards.
    # The loop ends when the parked pool is empty -- meter() returns 1 for that.
    sleep "$INTERVAL"
done
tail -5 "$OUT"
