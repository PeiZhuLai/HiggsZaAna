#!/usr/bin/env bash
# Resume chain_mhlow100.sh after step 4: the smoke test had passed (chunk0.root holds the limit tree), only the
# checker failed (PyROOT dropped the temporary TFile). Submit, watch, summarize, as the end of chain_mhlow100.sh.
# Output: bias_screen_ms_lo100_mA10.submit, out/ms_lo100/, resume_mhlow100.log (stdout), watch_ms_lo100.txt, summary_ms_lo100.txt
set -o pipefail
B=$(cd "$(dirname "$0")" && pwd); cd $B
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh); cd $B
python3 gen_submit.py ms_lo100 10 || { echo "STOP submit file"; exit 1; }
n=$(grep -c '^ms_lo100 ' bias_screen_ms_lo100_mA10.submit)
[ "$n" -le 120 ] || { echo "STOP $n jobs, more than the 120 announced"; exit 1; }
sub=$(env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=$KRB5CCNAME XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR condor_submit bias_screen_ms_lo100_mA10.submit)
echo "$sub"; cl=$(echo "$sub" | sed -n 's/.*submitted to cluster \([0-9]*\).*/\1/p'); [ -n "$cl" ] || { echo "STOP submit"; exit 1; }
echo "CLUSTER $cl jobs $n $(date +%F_%T)"
./watch.sh ms_lo100 $cl $n; rc=$?
[ $rc = 0 ] && python3 summarize.py ms_lo100 && cat summary_ms_lo100.txt
echo "END rc=$rc $(date +%F_%T)"; exit $rc
