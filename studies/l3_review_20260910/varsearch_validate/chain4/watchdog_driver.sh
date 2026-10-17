#!/usr/bin/env bash
# Detached watchdog of drive_step4.sh: if the driver process disappears without having written a DONE/STOP line,
# append DIED to the registry. Numbers only.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
REG=$V/logs/PENDING_RESULTS.txt
for _ in $(seq 1 2880); do
  al=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[d]rive_step4.sh' || true)
  if [ "$al" -eq 0 ]; then
    grep -q "step4" <(tail -3 $REG) || echo "$(date '+%F %H:%M') DIED  step4 driver  last status: $(cat $V/chain4/logs/drive_step4.status 2>/dev/null)" >> $REG
    exit 0
  fi
  sleep 60
done
