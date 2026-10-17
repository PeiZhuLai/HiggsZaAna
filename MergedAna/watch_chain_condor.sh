#!/usr/bin/env bash
# Detached watcher for one condor chain job: lands DONE/DIED in the registry.
# usage: watch_chain_condor.sh <cluster> <chain.log>
CL=$1; LOG=$2
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
for _ in $(seq 1 600); do
    st=$(condor_q "$CL" -af JobStatus 2>/dev/null | head -1)
    [ "$st" = 5 ] && { echo "$(date '+%F %H:%M') HELD  merged_chain $CL" >> "$REG"; exit 1; }
    if [ -z "$st" ]; then
        if grep -q 'chain exit code: 0' "$LOG" 2>/dev/null; then
            echo "$(date '+%F %H:%M') DONE  merged_chain $CL -> $LOG" >> "$REG"
        else
            echo "$(date '+%F %H:%M') DIED  merged_chain $CL last: $(grep '====' "$LOG" 2>/dev/null | tail -1)" >> "$REG"
        fi
        exit 0
    fi
    sleep 120
done
echo "$(date '+%F %H:%M') TIMEOUT merged_chain $CL" >> "$REG"
