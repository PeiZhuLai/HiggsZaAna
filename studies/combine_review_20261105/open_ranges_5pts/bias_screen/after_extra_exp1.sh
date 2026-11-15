#!/usr/bin/env bash
# Second batch of the mA10 bias screen (user-approved): wait for the fTest with EXTRA_ENV_PDFS=Exp1, then card,
# truths (every envelope member), submit file and condor_submit. Waits on the .done marker only, hard cap 3 h.
# Input : ../fit_results_ms_extra_Exp1/10/{.done,CMS-HGG_mva_13p6TeV_multipdf.root}
# Output: cards/ms_exp1/10_*, bias_screen_ms_exp1_mA10.submit, after_extra_exp1.log (cluster id), .after_exp1_done
set -o pipefail
B=$(cd "$(dirname "$0")" && pwd); S=$(dirname "$B"); R=$S/fit_results_ms_extra_Exp1/10
t0=$(date +%s)
while [ ! -f $R/.done ]; do
  [ $(( $(date +%s) - t0 )) -gt 10800 ] && { echo "TIMEOUT waiting for $R/.done"; exit 1; }
  sleep 60
done
echo "fTest: $(cat $R/.done)"; grep -q 'rc=0' $R/.done || { echo "fTest failed"; exit 1; }
grep 'EXTRA_ENV_PDFS' $R/ftest.stdout.log
grep -q 'EXTRA_ENV_PDFS: adding Exp1' $R/ftest.stdout.log || { echo "Exp1 was not added, not submitting"; exit 1; }
$B/make_card.sh 10 fit_results_ms_extra_Exp1 ms_exp1 || exit 1
( cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh) \
  && cd $B && python3 make_truths.py 10 ms_exp1 && python3 - <<'EOF'
import json; p='cards/ms_exp1/10_truths.json'; d=json.load(open(p)); d.pop('widePow1', None); json.dump(d, open(p,'w'), indent=1); print('truths', list(d))
EOF
) || exit 1
cd $B && python3 gen_submit.py ms_exp1 10 || exit 1
env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=$KRB5CCNAME XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR condor_submit bias_screen_ms_exp1_mA10.submit
date > $B/.after_exp1_done
