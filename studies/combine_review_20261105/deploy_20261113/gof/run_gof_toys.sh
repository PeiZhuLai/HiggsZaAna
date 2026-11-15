#!/bin/bash
# Sideband GOF toys (Baker-Cousins S, 300 toys, fit range from 95 GeV, SR bins dropped when reading) for every member
# of the redeployed envelopes of mA9, 11, 12, 13, 16, 17, 20, read from production fit_results_run3/<m> after stage 3.
# At most 5 processes. Script: ../../gof_fix_mA15_28/gof_toys.py (unchanged).
# Output: toys/mA<m>_<pdf>_lo95.json, logs/gof_<m>_<pdf>.log, completion marker logs/gof_toys.DONE
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$W/../../gof_fix_mA15_28/gof_toys.py
FR=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
rm -f logs/gof_toys.DONE; export G
for m in 9 11 12 13 16 17 20; do
  for p in $(grep -oE "pdf : [A-Za-z0-9]+" $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt | awk '{print $3}'); do echo "$m $p"; done
done | xargs -P 5 -L 1 bash -c 'python3 $G $0 $1 95 300 $((9000+$0*10)) > logs/gof_$0_$1.log 2>&1'
date > logs/gof_toys.DONE
