#!/usr/bin/env bash
# Registry watcher for the subset searches (systemd units hza-lowsubset-{backward,forward}). Numbers only.
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/lowparam
REG=$(dirname $W)/logs/PENDING_RESULTS.txt; left="backward forward"
for _ in $(seq 1 2880); do
  rest=""
  for M in $left; do
    if [ -e $W/state_subset_$M/DONE ]; then echo "$(date '+%F %H:%M') DONE  low subset-of-high search ($M)  -> $W/state_subset_$M/state.json" >> $REG
    elif ! systemctl --user is-active -q hza-lowsubset-$M; then echo "$(date '+%F %H:%M') DIED  low subset-of-high search ($M)  last: $(tail -1 $W/state_subset_$M/driver.log)" >> $REG
    else rest="$rest $M"; fi
  done
  left=$rest; [ -z "$left" ] && exit 0
  sleep 120
done
