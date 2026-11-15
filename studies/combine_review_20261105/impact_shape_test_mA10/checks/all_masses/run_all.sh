#!/bin/bash
# Expected limits with production and widened turn-on ranges at all 30 masses (limit_wide_ranges.py, --run blind).
# Five background processes, six masses each, every one in its own directory so that the combine files do not collide.
# Input : Combine/root_t2w_xsec/<m>_Datacard_leptons.root (read only)
# Output: chunk<i>/limit_wide_ranges.txt, chunk<i>/run.log, and the merged all_masses.txt written at the end
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
cd "$SRC" && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
i=0
for ms in "1 6 11 16 21 26" "2 7 12 17 22 27" "3 8 13 18 23 28" "4 9 14 19 24 29" "5 10 15 20 25 30"; do
  i=$((i+1)); mkdir -p "$W/chunk$i"; cd "$W/chunk$i" || exit 1
  python3 "$W/../limit_wide_ranges.py" $ms > run.log 2>&1 &
done
wait
cd "$W" || exit 1
{ head -1 chunk1/limit_wide_ranges.txt; cat chunk?/limit_wide_ranges.txt | grep -v '^mA' | sort -n; } > all_masses.txt
echo "[$(date +%H:%M:%S)] ALL DONE"
