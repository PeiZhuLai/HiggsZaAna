#!/bin/bash
# Waits for clusters 13420026 / 13420027 (mA1, mA2 bias screen with turn-on 95-113 GeV, 160 chunks).
# Completion read from out/cap113/*/chunk<c>.done only. Every 30 min; exits on all done, on any held job, or after 24 h.
# Status appended to logs/watch_status.txt; logs/ALL.DONE when complete.
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_mA1_2/bias
for i in $(seq 1 48); do
  a=$(find $R/out/cap113 -name 'chunk[0-9].done' 2>/dev/null | wc -l)
  q=$(condor_q 13420026 13420027 -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q 13420026 13420027 -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') cap113 $a/160  queue[$q] held $held" >> $R/logs/watch_status.txt
  [ "$a" -eq 160 ] && { echo ALL_DONE >> $R/logs/watch_status.txt; date > $R/logs/ALL.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $R/logs/watch_status.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $R/logs/watch_status.txt; exit 4
