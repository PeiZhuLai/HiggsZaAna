#!/bin/bash
# Waits for clusters 13422583 / 13422584 / 13422586 (mA11, 12 with 1 GeV width/smearing floor, mA28 with Bern5; 240 chunks).
# Completion read from out/flr/*/chunk<c>.done only. Every 30 min; exits on all done, on any held job, or after 24 h.
# Status appended to logs/watch_status.txt; logs/ALL.DONE when complete.
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_mA11_12_28/bias
for i in $(seq 1 48); do
  a=$(find $R/out/flr -name 'chunk[0-9].done' 2>/dev/null | wc -l)
  q=$(condor_q 13422583 13422584 13422586 -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q 13422583 13422584 13422586 -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') flr $a/240  queue[$q] held $held" >> $R/logs/watch_status.txt
  [ "$a" -eq 240 ] && { echo ALL_DONE >> $R/logs/watch_status.txt; date > $R/logs/ALL.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $R/logs/watch_status.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $R/logs/watch_status.txt; exit 4
