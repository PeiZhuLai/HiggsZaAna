#!/usr/bin/env bash
#
# Auto-release condor jobs held for exceeding their memory request, bumping the
# request one rung each time.
#
# WHY
#   The bkg + data re-production was submitted at the memory the original
#   production used (3000 MB). Measured on the live queue at 21:05 on 2026-09-12:
#       running jobs : median RSS 1657 MB, p90 5953 MB
#       held jobs    : n=766, mean 5635 MB, min 3672 MB, max 5980 MB
#   So the TYPICAL job is comfortable at 3000 MB and only a tail is not -- the
#   same shape the signal stage showed at fpo=1 (median 35 MB, max 4151 MB), just
#   4x up because bkg/data run at fpo=4. Raising RequestMemory for everything
#   would cost slots for the ~89 % of jobs that never need it, so instead each
#   job that actually hits the wall gets bumped individually and released.
#
#   This is not the fpo mistake repeating: fpo stays at 4. The knob being turned
#   here is memory for the jobs that measured over the limit, which is what the
#   measurement says to turn.
#
# WHY a loop and not a one-shot
#   Jobs keep reaching their memory peak as they run, so new holds appear for
#   hours. HZgamma's submit template sets PeriodicRelease = false, so nothing
#   releases them on its own and HiggsDNA retires a job after 5 failures --
#   a held job left alone is a lost chunk.
#
# LADDER
#   3000 -> 8000 -> 16000 -> 20000 MB. A job already at 20000 is left held and
#   reported: past that the slot essentially never matches and the right answer
#   is to rerun that file at fpo=1, not to keep bumping.
#
# USAGE
#   bash release_memory_holds.sh                      # dry run: report only
#   DRY_RUN=0 bash release_memory_holds.sh            # one pass
#   DRY_RUN=0 LOOP=1 nohup bash release_memory_holds.sh >/dev/null 2>&1 &
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
LOOP="${LOOP:-0}"
INTERVAL="${INTERVAL:-600}"
OUT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/memory_holds.log"
# Cmd-based, never batch name: another project's DY jobs share this schedd.
CONST='Owner=="pelai" && regexp("HZa/HiggsZaAna",Cmd)'

next_rung() {
    case "$1" in
        ''|*[!0-9]*) echo 8000 ;;
        *) if   [ "$1" -lt  8000 ]; then echo 8000
           elif [ "$1" -lt 16000 ]; then echo 16000
           elif [ "$1" -lt 20000 ]; then echo 20000
           else echo stuck; fi ;;
    esac
}

pass() {
    local ts; ts=$(date '+%F %T')
    # HoldReasonCode 34 == "job exceeded its memory request". Other hold codes
    # (input errors, signals) must NOT be released blindly -- they would just
    # burn another of the job's five allowed failures.
    local held; held=$(condor_q -constraint "$CONST && JobStatus==5" -af ClusterId ProcId HoldReasonCode RequestMemory 2>/dev/null \
        | awk '$3 == 34 {print $1"."$2, $4}')
    if [ -z "$held" ]; then
        echo "$ts  no memory holds" >> "$OUT"; return 0
    fi
    local n=0 stuck=0
    declare -A bump=()
    while read -r id mem; do
        local nm; nm=$(next_rung "$mem")
        if [ "$nm" = stuck ]; then stuck=$((stuck+1)); continue; fi
        bump[$nm]="${bump[$nm]:-} $id"
        n=$((n+1))
    done <<< "$held"
    echo "$ts  memory holds: ${n} to bump, ${stuck} already at 20000 (left held)" >> "$OUT"
    for nm in "${!bump[@]}"; do
        local ids=${bump[$nm]}
        local cnt; cnt=$(wc -w <<< "$ids")
        echo "$ts    -> RequestMemory=${nm} for ${cnt} jobs" >> "$OUT"
        if [ "$DRY_RUN" = "1" ]; then
            echo "    would run: condor_qedit <${cnt} ids> RequestMemory ${nm} && condor_release <same>"
        else
            # condor_qedit takes ids in batches; keep them small so one bad id
            # does not abort the whole set.
            xargs -n 40 <<< "$ids" | while read -r batch; do
                condor_qedit $batch RequestMemory "${nm}" >> "$OUT" 2>&1
                condor_release $batch >> "$OUT" 2>&1
            done
        fi
    done
}

while :; do
    pass
    [ "$LOOP" = "1" ] && [ "$DRY_RUN" != "1" ] || break
    # Stop once no driver is left: nothing will be resubmitted anyway.
    # Not tied to a driver being alive: the production ran driverless for most of
    # 2026-09-13 and this loop exited at 15:51, leaving 214 jobs stranded on
    # memory holds all night. Exit when there is nothing left in the queue.
    # An empty queue alone is not a reason to stop: a driver that is still
    # writing 28k job configs has nothing queued yet, and the releaser exiting
    # there leaves the whole production unattended (it did exactly that at
    # 22:25). Stop only when the queue is empty AND no driver is left to fill it.
    if [ "$(condor_q -constraint "$CONST" -af ClusterId 2>/dev/null | wc -l)" -eq 0 ] \
       && ! ps -u "$(id -un)" -o cmd= 2>/dev/null | awk '/run_analysis\.py/ && !/awk/' | grep -q .; then
        echo "$(date '+%F %T')  queue empty and no drivers -- releaser exiting" >> "$OUT"; break
    fi
    sleep "$INTERVAL"
done
tail -6 "$OUT"
