#!/usr/bin/env bash
# mA10 with the lower fit edge at 100 GeV instead of 95 (user-approved test): the fitted turn-on sits at 96-98 GeV,
# 1-3 GeV from the 95 GeV edge, with the width on its 0.1 floor, so the turn-on is not constrained by the data;
# mA1 already uses 100-180 GeV. Same multi-start binary and opened ranges as fit_results_ms, no EXTRA/EXCLUDE.
#  1 fTest MHLOW=100 (bin/fTest_ALP_turnOn.ms)
#  2 checks: rc, dataset reduced to [100,180], envelope, card observable 80 bins
#  3 card + truths (every envelope member at its stored values; injection r = 0 and the same q50 0.0199)
#  4 local smoke test (5 toys, first truth, tag ms_lo100_smoke), then condor <members> x 2 inj x 10 chunks x 1000 toys
#  5 watch.sh until all chunks are done or a job is held, then summarize.py
# Input : /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M10/ws/run3.root
# Output: ../fit_results_ms_lo100/10/, cards/ms_lo100/, out/ms_lo100/, out/ms_lo100_smoke/,
#         bias_screen_ms_lo100_mA10.submit, chain_mhlow100.log (stdout), watch_ms_lo100.txt, summary_ms_lo100.txt
set -o pipefail
B=$(cd "$(dirname "$0")" && pwd); S=$(dirname "$B"); R=$S/fit_results_ms_lo100/10
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
echo "== 1 fTest $(date +%F_%T)"
MS=1 MHLOW=100 RES_OVERRIDE=fit_results_ms_lo100 $S/run_ftest_v3.sh 10
echo "fTest: $(cat $R/.done) $(date +%F_%T)"; grep -q 'rc=0' $R/.done || { echo "STOP fTest failed"; exit 1; }
echo "== 2 checks"
grep -m1 'dataset reduced' $R/ftest.stdout.log || { echo "STOP dataset not reduced to the 100 GeV edge"; exit 1; }
grep 'multi-start' $R/ftest.stdout.log
grep 'pdf :' $R/HZAmassInde_fTest/EnvelopeResults.txt
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
$B/make_card.sh 10 fit_results_ms_lo100 ms_lo100 || { echo "STOP card"; exit 1; }
cd $B && python3 make_truths.py 10 ms_lo100 || { echo "STOP truths"; exit 1; }
python3 - <<'EOF' || { echo "STOP truth/card check"; exit 1; }
import json, ROOT
ROOT.gROOT.SetBatch(True)
p = 'cards/ms_lo100/10_truths.json'; d = json.load(open(p)); d = {k: v for k, v in d.items() if not k.startswith('wide')}
json.dump(d, open(p, 'w'), indent=1); print('truths', list(d))
w = ROOT.TFile.Open('cards/ms_lo100/10_Datacard_leptons.root').Get('w'); m = w.var('CMS_hza_mass')
print('card observable', m.getMin(), m.getMax(), m.getBins())
data = w.data('data_obs'); print('data_obs entries', data.numEntries(), 'sum', data.sumEntries())
assert abs(m.getMin() - 100.) < 1e-6 and m.getBins() == 80, 'observable is not 100-180 GeV in 80 bins'
EOF
first=$(python3 -c "import json; print(list(json.load(open('cards/ms_lo100/10_truths.json')))[0])")
echo "== 4 smoke $(date +%F_%T) truth $first"
ln -sfn ms_lo100 cards/ms_lo100_smoke
./runjob.sh ms_lo100_smoke 10 $first 0 0 5 || { echo "STOP smoke failed"; exit 1; }
python3 - <<EOF || { echo "STOP smoke content"; exit 1; }
import ROOT
f = ROOT.TFile.Open('out/ms_lo100_smoke/mA10_${first}_0/chunk0.root'); t = f.Get('limit'); print('smoke entries', t.GetEntries()); assert t.GetEntries() >= 3
EOF
python3 gen_submit.py ms_lo100 10 || { echo "STOP submit file"; exit 1; }
n=$(grep -c '^ms_lo100 ' bias_screen_ms_lo100_mA10.submit)
[ "$n" -le 120 ] || { echo "STOP $n jobs, more than the 120 announced"; exit 1; }
sub=$(env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=$KRB5CCNAME XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR condor_submit bias_screen_ms_lo100_mA10.submit)
echo "$sub"; cl=$(echo "$sub" | sed -n 's/.*submitted to cluster \([0-9]*\).*/\1/p'); [ -n "$cl" ] || { echo "STOP submit"; exit 1; }
echo "CLUSTER $cl jobs $n $(date +%F_%T)"
echo "== 5 watch"
./watch.sh ms_lo100 $cl $n; rc=$?
[ $rc = 0 ] && python3 summarize.py ms_lo100 && cat summary_ms_lo100.txt
echo "END rc=$rc $(date +%F_%T)"; exit $rc
