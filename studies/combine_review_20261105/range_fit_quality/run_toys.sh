#!/bin/bash
# toys_dS.py for the selected function at five masses, one background process each, 100 toys each.
# Output: toys/mA<m>_<pdf>.json, logs/toys_mA<m>.log
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
cd "$SRC" && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
for x in "10 Pow1" "15 Pow1" "20 Pow1" "22 Exp3" "23 Exp3"; do
  set -- $x; python3 toys_dS.py "$1" "$2" 100 $((1000+$1)) > logs/toys_mA$1.log 2>&1 &
done
wait; echo "[$(date +%H:%M:%S)] ALL DONE"
