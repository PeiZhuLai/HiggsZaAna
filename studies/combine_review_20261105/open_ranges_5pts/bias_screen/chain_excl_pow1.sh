#!/usr/bin/env bash
# mA10 without Pow1 (user decision: Pow1 truth fails the bias criterion, remove it from the envelope).
#  1 fTest with the multi-start binary + EXCLUDE_ENV_PDFS=Pow1 (bin/fTest_ALP_turnOn.msxe)
#  2 checks: Pow1 removed, envelope = Lau1, Exp3, Bern5, their stored values identical to fit_results_ms (same fits)
#  3 card + truths: the 3 members, plus Pow1 at its fit_results_ms values generated from cards/ms (gen_card),
#    to show what a Pow1-like background does to the envelope without it
#  4 local smoke test (5 toys, truth Pow1, separate tag ms_noPow1_smoke), then condor 4 truths x 2 inj x 10 = 80 jobs
#  5 watch.sh until all chunks are done or a job is held
# Seeds follow runjob.sh, so every truth uses the same toys as the ms batch (paired comparison).
# Input : /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M10/ws/run3.root, cards/ms/10_{truths.json,Datacard_leptons.root}
# Output: ../fit_results_ms_excl_Pow1/10/, cards/ms_noPow1/, out/ms_noPow1/, out/ms_noPow1_smoke/,
#         bias_screen_ms_noPow1_mA10.submit, chain_excl_pow1.log (this script's stdout), watch_ms_noPow1.txt
set -o pipefail
B=$(cd "$(dirname "$0")" && pwd); S=$(dirname "$B"); R=$S/fit_results_ms_excl_Pow1/10
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
echo "== 1 fTest $(date +%F_%T)"
MS=1 BIN_OVERRIDE=./bin/fTest_ALP_turnOn.msxe RES_OVERRIDE=fit_results_ms_excl_Pow1 EXCLUDE_ENV_PDFS=Pow1 $S/run_ftest_v2.sh 10
echo "fTest: $(cat $R/.done) $(date +%F_%T)"; grep -q 'rc=0' $R/.done || { echo "STOP fTest failed"; exit 1; }
echo "== 2 checks"
grep 'EXCLUDE_ENV_PDFS' $R/ftest.stdout.log
grep -q 'EXCLUDE_ENV_PDFS: removing Pow1' $R/ftest.stdout.log || { echo "STOP Pow1 not removed"; exit 1; }
grep 'pdf :' $R/HZAmassInde_fTest/EnvelopeResults.txt
cd $FF/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
$B/make_card.sh 10 fit_results_ms_excl_Pow1 ms_noPow1 || { echo "STOP card"; exit 1; }
cd $B && python3 make_truths.py 10 ms_noPow1 || { echo "STOP truths"; exit 1; }
python3 - <<'EOF' || { echo "STOP truth check"; exit 1; }
import json
new = json.load(open('cards/ms_noPow1/10_truths.json')); old = json.load(open('cards/ms/10_truths.json'))
new.pop('widePow1', None)
assert sorted(new) == ['Bern5', 'Exp3', 'Lau1'], 'envelope is %s' % sorted(new)
for k in new: assert new[k]['set'] == old[k]['set'], '%s values differ: %s vs %s' % (k, new[k]['set'], old[k]['set'])
p = dict(old['Pow1']); p['gen_card'] = '/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/open_ranges_5pts/bias_screen/cards/ms/10_Datacard_leptons.root'
new['Pow1'] = p
json.dump(new, open('cards/ms_noPow1/10_truths.json', 'w'), indent=1)
print('truths', {k: (v['index'], v.get('gen_card', 'own card')) for k, v in new.items()}); print('member values identical to fit_results_ms')
EOF
echo "== 4 smoke $(date +%F_%T)"
ln -sfn ms_noPow1 cards/ms_noPow1_smoke
./runjob_v2.sh ms_noPow1_smoke 10 Pow1 0 0 5 || { echo "STOP smoke failed"; exit 1; }
[ -f out/ms_noPow1_smoke/mA10_Pow1_0/chunk0.done ] || { echo "STOP smoke no output"; exit 1; }
python3 - <<'EOF' || { echo "STOP smoke content"; exit 1; }
import ROOT
f = ROOT.TFile.Open('out/ms_noPow1_smoke/mA10_Pow1_0/chunk0.root'); t = f.Get('limit'); n = t.GetEntries()
print('smoke entries', n); assert n >= 3
EOF
EXE=runjob_v2.sh python3 gen_submit.py ms_noPow1 10 || { echo "STOP submit file"; exit 1; }
sub=$(env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=$KRB5CCNAME XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR condor_submit bias_screen_ms_noPow1_mA10.submit)
echo "$sub"; cl=$(echo "$sub" | sed -n 's/.*submitted to cluster \([0-9]*\).*/\1/p'); [ -n "$cl" ] || { echo "STOP submit"; exit 1; }
echo "CLUSTER $cl $(date +%F_%T)"
echo "== 5 watch"
./watch.sh ms_noPow1 $cl 80; rc=$?
[ $rc = 0 ] && python3 summarize.py ms_noPow1 && cat summary_ms_noPow1.txt
echo "END rc=$rc $(date +%F_%T)"; exit $rc
