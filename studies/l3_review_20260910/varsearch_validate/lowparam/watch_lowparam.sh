#!/usr/bin/env bash
# Registry watcher for the two low-only searches (systemd units hza-lowparam-{shared,free}): DONE when state_<p>/DONE
# exists, DIED when the unit is no longer active without it. Numbers only.
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/lowparam
REG=$(dirname $W)/logs/PENDING_RESULTS.txt; left="shared free"
for _ in $(seq 1 2880); do
  rest=""
  for P in $left; do
    if [ -e $W/state_$P/DONE ]; then echo "$(date '+%F %H:%M') DONE  lowparam search ($P pool)  -> $W/state_$P/state.json" >> $REG
    elif ! systemctl --user is-active -q hza-lowparam-$P; then echo "$(date '+%F %H:%M') DIED  lowparam search ($P pool)  last: $(tail -1 $W/state_$P/driver.log)" >> $REG
    else rest="$rest $P"; fi
  done
  left=$rest; [ -z "$left" ] && exit 0
  sleep 120
done
