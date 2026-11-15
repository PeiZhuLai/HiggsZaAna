#!/bin/bash
# Re-runs the members of logs/members.txt whose toys/<tag>_an.json is missing (the first pass of run_gof_toys.sh lost
# mA17-30 when the AFS work quota was full: the log redirection failed, so the command never started). Same settings.
# Output: toys/mA<m>_<pdf>_lo<lo>_an.json, logs/gof_<m>_<pdf>.log, logs/rc_missing.txt, logs/gof_missing.DONE
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$W/../gof_fix_mA15_28/gof_toys.py
FR=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
cd "$W" || exit 1
rm -f logs/gof_missing.DONE; export G FR
while read m p lo; do [ -s toys/mA${m}_${p}_lo${lo}_an.json ] || echo "$m $p $lo"; done < logs/members.txt > logs/members_missing.txt
xargs -P 4 -L 1 bash -c 'MULTIPDF=$FR/$0/CMS-HGG_mva_13p6TeV_multipdf.root TAG=_an python3 $G $0 $1 $2 300 $((9700+$0*10)) > logs/gof_$0_$1.log 2>&1; echo "$0 $1 rc=$?" >> logs/rc_missing.txt' < logs/members_missing.txt
date > logs/gof_missing.DONE
