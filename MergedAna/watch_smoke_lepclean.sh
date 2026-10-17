#!/usr/bin/env bash
# Detached watcher for the lepclean data-friend smoke test driver.
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/smoke_lepclean_Muon1_RunD.log
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
for _ in $(seq 1 360); do
    if grep -q '^Done: ' "$LOG" 2>/dev/null; then
        echo "$(date '+%F %H:%M') DONE  smoke_lepclean_Muon1_RunD  $(grep '^Done: ' "$LOG")" >> "$REG"; exit 0
    fi
    alive=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_merged_data_ml_friend.sh' || true)
    if [ "$alive" -eq 0 ]; then
        echo "$(date '+%F %H:%M') DIED  smoke_lepclean_Muon1_RunD driver gone without Done line" >> "$REG"; exit 1
    fi
    sleep 120
done
echo "$(date '+%F %H:%M') TIMEOUT smoke_lepclean_Muon1_RunD" >> "$REG"
