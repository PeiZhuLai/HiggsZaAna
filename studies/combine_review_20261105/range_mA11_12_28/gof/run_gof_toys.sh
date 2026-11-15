#!/bin/bash
# Sideband GOF toys (300 toys, SR bins dropped when reading) for every member of the staging envelopes of mA11, 12
# (width/smearing floor 1 GeV) and mA28 (Bernstein up to order 6), fit from 95 GeV. Script
# ../../gof_fix_mA15_28/gof_toys.py with MULTIPDF=<staging multipdf>, TAG=_flr.
# Output: toys/mA<m>_<pdf>_lo95_flr.json, logs/gof_<m>_<pdf>.log, logs/gof_toys.DONE
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$W/../../gof_fix_mA15_28/gof_toys.py; R=$W/../fit_results
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
rm -f logs/gof_toys.DONE; export G R
for m in 11 12 28; do
  for p in $(grep -oE "pdf : [A-Za-z0-9]+" $R/$m/HZAmassInde_fTest/EnvelopeResults.txt | awk '{print $3}'); do echo "$m $p 95"; done
done | xargs -P 5 -L 1 bash -c 'MULTIPDF=$R/$0/CMS-HGG_mva_13p6TeV_multipdf.root TAG=_flr python3 $G $0 $1 $2 300 $((9700+$0*10)) > logs/gof_$0_$1.log 2>&1'
date > logs/gof_toys.DONE
