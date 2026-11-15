#!/bin/bash
# Sideband GOF toys: mA28 (Exp1, Pow1, Lau1) x lower edge 95/100/105 GeV, mA15 (Pow1, Lau1, Exp1, Bern5) at 95 GeV.
# 300 toys each, at most 6 processes at a time. Output: toys/*.json, logs/gof_*.log, completion marker logs/gof_toys.DONE
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
rm -f logs/gof_toys.DONE
{ for lo in 95 100 105; do for p in Exp1 Pow1 Lau1; do echo "28 $p $lo"; done; done
  for p in Pow1 Lau1 Exp1 Bern5; do echo "15 $p 95"; done; } |
xargs -P 6 -L 1 bash -c 'python3 gof_toys.py $0 $1 $2 300 $((7000+$0*10+${2})) > logs/gof_$0_$1_lo$2.log 2>&1'
date > logs/gof_toys.DONE
