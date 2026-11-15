#!/bin/bash
# fit_quality.py at all 30 masses, five background processes. Sidebands only.
# Input : Background/ALP_BkgModel_ReReco/fit_results_run3/<mA>/CMS-HGG_mva_13p6TeV_multipdf.root (read only)
# Output: results/mA<mA>.json, logs/chunk<i>.log
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
cd "$SRC" && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" && mkdir -p logs results; i=0
for ms in "1 6 11 16 21 26" "2 7 12 17 22 27" "3 8 13 18 23 28" "4 9 14 19 24 29" "5 10 15 20 25 30"; do
  i=$((i+1)); ( for m in $ms; do python3 fit_quality.py "$m" results; done ) > logs/chunk$i.log 2>&1 &
done
wait; echo "[$(date +%H:%M:%S)] ALL DONE"
