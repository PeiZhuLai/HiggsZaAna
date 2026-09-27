#!/usr/bin/env bash
#
# Persistent watchdog for the B-mode bkg + data re-production.
#
# WHY a script and not an interactive check:
#   The two hangs that cost the most time today -- the schedd refusing queries
#   (SECMAN:2011, 17:49-20:51) and the driver stuck in do_wait with a frozen log
#   -- were BOTH invisible to "is the process alive" and to "did the driver
#   return 0" (run_analysis.py counts retired jobs as complete and returns 0 even
#   when every job died). The only signal that caught them was the chunk count
#   not moving, so that delta is the primary metric here.
#
# WHY there is no grep over eos_logs:
#   eos_logs holds hundreds of thousands of files from every past production --
#   `ls | head` on it does not return within two minutes over FUSE. A watchdog
#   that walks it takes longer than its own interval and starves itself. Held
#   jobs and their hold reasons come from condor_q instead, which is O(queue).
#
# USAGE
#   bash watch_bkgdata.sh              # dry run: one cycle, prints, exits
#   DRY_RUN=0 nohup bash watch_bkgdata.sh >/dev/null 2>&1 &
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
INTERVAL="${INTERVAL:-1800}"
D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
S="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix"
OUT="${D}/logs_fsrfix/watch_bkgdata.log"

count() { find "$1" -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l | tr -d ' '; }

pb=-1; pd=-1
while :; do
    ts=$(date '+%F %T')
    nb=$(count "$S/Bkg_MC"); nd=$(count "$S/Data")
    # ps|awk, never pgrep -f: pgrep would match this script's own command line.
    drv=$(ps -u "$(id -un)" -o cmd= 2>/dev/null | awk '/run_analysis\.py/ && !/awk/' | wc -l)
    # Parse by LABEL, not by field position. The totals line is
    #   "Total for pelai: N jobs; N completed, N removed, N idle, N running, N held, N suspended"
    # so positional awk silently hands back the words "completed,"/"removed,"
    # instead of numbers, and every numeric test then compares strings.
    read -r tot idle run held < <(condor_q -totals 2>/dev/null | sed -n \
        's/.*Total for pelai: \([0-9]*\) jobs; [0-9]* completed, [0-9]* removed, \([0-9]*\) idle, \([0-9]*\) running, \([0-9]*\) held.*/\1 \2 \3 \4/p')
    left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
    if [ "$pb" -lt 0 ]; then db=init; dd=init; else db=$((nb-pb)); dd=$((nd-pd)); fi
    printf '%s bkg=%s(+%s) data=%s(+%s) drivers=%s queue=%s idle=%s run=%s held=%s proxy=%ss\n' \
        "$ts" "$nb" "$db" "$nd" "$dd" "$drv" "${tot:-0}" "${idle:-0}" "${run:-0}" "${held:-0}" "$left" >> "$OUT"

    [ "${left:-0}" -lt 3600 ] && echo "$ts  ALERT proxy under 1 h (${left}s) -- jobs starting after expiry lose xrootd" >> "$OUT"
    if [ "${held:-0}" -gt 0 ] 2>/dev/null; then
        echo "$ts  ALERT ${held} held -- top hold reasons:" >> "$OUT"
        condor_q -held -af HoldReason 2>/dev/null | cut -c1-110 | sort | uniq -c | sort -rn | head -3 \
            | sed "s/^/$ts    /" >> "$OUT"
    fi
    if [ "$db" != init ] && [ "$db" -eq 0 ] && [ "$dd" -eq 0 ] && [ "${run:-0}" -gt 0 ]; then
        echo "$ts  ALERT ${run} jobs running but no new chunk in ${INTERVAL}s -- suspect a stalled driver or schedd" >> "$OUT"
    fi
    # Do NOT exit when the drivers are gone. From 2026-09-13 15:45 the production
    # runs deliberately driverless: HiggsDNA's retry path calls condor_release,
    # which drained the throttle's parked pool and put 2952 jobs back in flight.
    # The queue itself is the thing to watch now; it ends when the queue empties.
    [ "${tot:-0}" -eq 0 ] && { echo "$ts  queue empty -- watchdog exiting" >> "$OUT"; break; }

    pb=$nb; pd=$nd
    [ "$DRY_RUN" = "1" ] && { echo "(dry run: one cycle only)"; tail -2 "$OUT"; break; }
    sleep "$INTERVAL"
done
