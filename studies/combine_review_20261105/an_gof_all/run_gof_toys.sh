#!/bin/bash
# Sideband goodness of fit with toys for every member of the production envelope at the 30 mass points, for the AN
# table: 300 toys per member, SR bins (115-135 GeV) dropped when the data are read, fit from 95 GeV (100 GeV at mA1,
# as in production). Script ../gof_fix_mA15_28/gof_toys.py on the production multipdf (read-only).
# Input : flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/<m>/{CMS-HGG_mva_13p6TeV_multipdf.root,HZAmassInde_fTest/EnvelopeResults.txt}
# Output: toys/mA<m>_<pdf>_lo<lo>_an.json, logs/gof_<m>_<pdf>.log, logs/gof_toys.DONE
# At most 4 local processes.
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$W/../gof_fix_mA15_28/gof_toys.py
FR=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
mkdir -p logs toys; rm -f logs/gof_toys.DONE; export G FR
for m in $(seq 1 30); do
  lo=95; [ $m = 1 ] && lo=100
  for p in $(grep -oE "pdf : [A-Za-z0-9]+" $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt | awk '{print $3}'); do echo "$m $p $lo"; done
done > logs/members.txt
xargs -P 4 -L 1 bash -c 'MULTIPDF=$FR/$0/CMS-HGG_mva_13p6TeV_multipdf.root TAG=_an python3 $G $0 $1 $2 300 $((9700+$0*10)) > logs/gof_$0_$1.log 2>&1; echo "$0 $1 rc=$?" >> logs/rc.txt' < logs/members.txt
date > logs/gof_toys.DONE
