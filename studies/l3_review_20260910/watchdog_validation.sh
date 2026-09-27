#!/usr/bin/env bash
#
# Detached watchdog for the copy-then-read (staging) validation batch.
#
# WHY DETACHED + A REGISTRY
#   Local background processes in this environment get killed without notice
#   (six times in this session). A detached watchdog lands the result on disk,
#   and the registry line is what makes it reportable: the result must not sit
#   on disk waiting to be asked for.
#
# It writes NUMBERS ONLY. The verdict is not written here -- the same counts can
# mean different things (e.g. exit 127 in 3 s is a broken submit environment, not
# a failed code change, which is exactly what happened at 12:39 today).
#
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
S=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix
SMP=DYJetsToLL_2022postEE
REG=$D/logs_fsrfix/PENDING_RESULTS.txt
STATUS=$D/logs_fsrfix/validation_status.txt
RES=$D/logs_fsrfix/validation_result.txt
CONST='Owner=="pelai" && regexp("HZa/HiggsZaAna",Cmd)'
BASE=235          # chunks in this sample before the batch was submitted
CUT=$(date -d '2026-09-15 17:24' +%s)

for i in $(seq 1 400); do
    n=$(find $S/Bkg_MC/$SMP -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l)
    q=$(condor_q -constraint "$CONST" -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
    echo "poll $(date '+%F %H:%M:%S') chunks=$n (+$((n-BASE))) queue=[$q]" > "$STATUS"
    left=$(condor_q -constraint "$CONST" -af ClusterId 2>/dev/null | wc -l)
    if [ "$left" -eq 0 ]; then
        {
          echo "=== staging validation, $SMP, submitted 17:24 $(date '+%F') ==="
          echo "chunks before=$BASE after=$n delta=$((n-BASE)) of 20 submitted"
          condor_history -limit 120 -const "$CONST" -af CompletionDate ExitCode RemoteWallClockTime Err RemoveReason 2>/dev/null \
            | awk -v cut=$CUT '$1+0>cut {r=($0~/wall time/)?"WALL":"exit"$2; c[r]++}
                 END{for(k in c) printf "outcome %-10s n=%d\n", k, c[k]}'
          echo "--- error classes among failures ---"
          condor_history -limit 120 -const "$CONST" -af CompletionDate ExitCode Err 2>/dev/null \
            | awk -v cut=$CUT '$1+0>cut && $2!=0 {print $3}' | while read -r e; do
                [ -f "$e" ] || continue
                if   grep -q 'vector_read properly' "$e" 2>/dev/null; then echo "vector_read"
                elif grep -q -e '\[3010\]' -e 'Unable to open' "$e" 2>/dev/null; then echo "unavailable"
                elif grep -q 'xrdcp failed' "$e" 2>/dev/null; then echo "xrdcp_failed"
                else echo "other"; fi
              done | sort | uniq -c
        } > "$RES" 2>&1
        echo "$(date '+%F %H:%M') DONE  staging-validation  -> $RES" >> "$REG"
        exit 0
    fi
    sleep 120
done
echo "$(date '+%F %H:%M') TIMEOUT  staging-validation  after 400 polls" >> "$REG"
