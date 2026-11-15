#!/bin/bash
# toys_dS.py for mA28 Exp1 (worst sideband S of the 30 masses; ranges already wide, so only S_prod matters), 100 toys.
# Output: toys/mA28_Exp1.json, logs/toys_mA28.log
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
cd "$SRC" && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
python3 toys_dS.py 28 Exp1 100 1028 > logs/toys_mA28.log 2>&1
echo "[$(date +%H:%M:%S)] ALL DONE"
