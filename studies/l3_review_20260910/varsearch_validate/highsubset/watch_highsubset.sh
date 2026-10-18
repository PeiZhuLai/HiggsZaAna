#!/usr/bin/env bash
# Registry watcher for the high-mass window-training subset search (systemd unit hza-highsubset-backward).
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/highsubset
REG=$(dirname $H)/logs/PENDING_RESULTS.txt
echo "$(date '+%F %H:%M') START high subset (window training) backward search -> $H/state_backward/state.json" >> $REG
for _ in $(seq 1 2880); do
  if [ -e $H/state_backward/DONE ]; then echo "$(date '+%F %H:%M') DONE  high subset backward -> $H/state_backward/state.json" >> $REG; exit 0; fi
  systemctl --user is-active -q hza-highsubset-backward || { echo "$(date '+%F %H:%M') DIED  high subset backward: $(journalctl --user -u hza-highsubset-backward --no-pager | tail -1)" >> $REG; exit 1; }
  sleep 120
done
