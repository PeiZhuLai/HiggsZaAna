#!/bin/bash
# Sideband GOF toys (300 toys, SR bins dropped when reading) for every member of the mA11 envelope without the 1 GeV
# floor, fit from 95 GeV. Script ../../gof_fix_mA15_28/gof_toys.py with MULTIPDF=<staging multipdf>, TAG=_nf.
# Output: toys/mA11_<pdf>_lo95_nf.json, logs/gof_11_<pdf>.log, logs/gof_toys.DONE
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$W/../../gof_fix_mA15_28/gof_toys.py; R=$W/../fit_results
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
mkdir -p logs toys; rm -f logs/gof_toys.DONE; export G R
for p in $(grep -oE "pdf : [A-Za-z0-9]+" $R/11/HZAmassInde_fTest/EnvelopeResults.txt | awk '{print $3}'); do echo "11 $p 95"; done \
 | xargs -P 4 -L 1 bash -c 'MULTIPDF=$R/$0/CMS-HGG_mva_13p6TeV_multipdf.root TAG=_nf python3 $G $0 $1 $2 300 $((9700+$0*10)) > logs/gof_$0_$1.log 2>&1'
date > logs/gof_toys.DONE
