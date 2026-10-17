#!/usr/bin/env bash
# Registry watcher for chain4.sh (restarted at stage 3 on 2026-10-14 after the AFS-quota incident): DONE / FAILED / DIED
# line in ../logs/PENDING_RESULTS.txt; also stops chain4 if any of its logs exceeds 200 MB (runaway output). Numbers only.
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chain4
REG=$(dirname $CH)/logs/PENDING_RESULTS.txt
for _ in $(seq 1 4320); do
  if [ -e $CH/logs/chain4.DONE ]; then echo "$(date '+%F %H:%M') DONE  chain4 (stages 3-10)  -> $CH/logs/chain4.log" >> $REG; exit 0; fi
  if [ -e $CH/logs/chain4.FAILED ]; then echo "$(date '+%F %H:%M') STOP  chain4  $(cat $CH/logs/chain4.FAILED)" >> $REG; exit 1; fi
  big=$(find $CH/logs $CH/flashggFinalFit -name '*.log' -size +200M 2>/dev/null | head -1)
  if [ -n "$big" ]; then
    for p in $(ps -u "$USER" -o pid=,cmd= | awk '/[c]hain4\.sh/ {print $1}'); do kill $p; done
    echo "$(date '+%F %H:%M') STOP  chain4 killed: runaway log $big ($(stat -c %s $big) bytes)" >> $REG; exit 1
  fi
  if [ "$(ps -u "$USER" -o cmd= | grep -c '[c]hain4\.sh')" -eq 0 ]; then echo "$(date '+%F %H:%M') DIED  chain4  last: $(tail -1 $CH/logs/chain4_heartbeat.txt)" >> $REG; exit 1; fi
  sleep 60
done
