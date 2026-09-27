#!/usr/bin/env bash
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while true; do
  act=$(systemctl --user is-active hza-train 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) unit=$act"
    echo "  stage    : $(grep -o '\[[0-9]/3\][^=]*' $D/logs_fsrfix/train.log 2>/dev/null | tail -1)"
    echo "  log lines: $(wc -l < $D/logs_fsrfix/train.log 2>/dev/null)"
    echo "  last     : $(tail -1 $D/logs_fsrfix/train.log 2>/dev/null | cut -c1-130)"
  } > $D/logs_fsrfix/train_heartbeat.txt
  [ "$act" != "active" ] && { echo "  UNIT FINISHED" >> $D/logs_fsrfix/train_heartbeat.txt; break; }
  sleep 600
done
