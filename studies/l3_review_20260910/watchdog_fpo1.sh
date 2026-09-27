#!/usr/bin/env bash
# Detached watchdog for the fpo=1 re-production. Numbers only, no verdicts.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
NEW=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1
REG=$D/logs_fsrfix/PENDING_RESULTS.txt
ST=$D/logs_fsrfix/fpo1_status.txt
CONST='Owner=="pelai" && regexp("HZa/HiggsZaAna",Cmd)'
TARGET=27898
prev=-1; stall=0
while :; do
    n=$(find $NEW -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l)
    q=$(condor_q -constraint "$CONST" -af JobStatus 2>/dev/null | sort | uniq -c | tr -d '\n' | tr -s ' ')
    cpu=$(condor_q -constraint "$CONST" -af RequestCpus 2>/dev/null | sort | uniq -c | tr -d '\n' | tr -s ' ')
    # Match on the OUTPUT DIR, not just run_analysis.py: other projects run the
    # same script on this login node (an HZgamma production showed up as two
    # extra "drivers" at 00:14 on 2026-09-16). Same class of mistake as an
    # unanchored batch-name regex.
    drv=$(ps -u "$(id -un)" -o cmd= 2>/dev/null | awk '/run_analysis\.py/ && /parquet_DNA_tmp_fsrfix_fpo1/ && !/awk/' | wc -l)
    echo "$(date '+%F %H:%M:%S') chunks=$n/$TARGET queue=[$q] cpus=[$cpu] drivers=$drv" > "$ST"
    if [ "$prev" -ge 0 ] && [ "$n" -eq "$prev" ]; then stall=$((stall+1)); else stall=0; fi
    prev=$n
    if [ "$n" -ge "$TARGET" ]; then
        echo "$(date '+%F %H:%M') DONE  fpo1-production  chunks=$n/$TARGET" >> "$REG"; exit 0; fi
    if [ "$drv" -eq 0 ] && [ "$(condor_q -constraint "$CONST" -af ClusterId 2>/dev/null | wc -l)" -eq 0 ]; then
        echo "$(date '+%F %H:%M') STOPPED  fpo1-production  chunks=$n/$TARGET drivers=0 queue=0" >> "$REG"; exit 1; fi
    [ "$stall" -ge 12 ] && { echo "$(date '+%F %H:%M') STALLED  fpo1-production  chunks=$n for 2h" >> "$REG"; stall=0; }
    sleep 600
done
