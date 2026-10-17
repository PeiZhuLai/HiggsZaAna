#!/usr/bin/env bash
# Re-run 10-17 14:5x: the first rescore failed the gate of rescore_selected.py (old low-mass features vs the deployed relE2 model); gate fixed.
# Full m_a 1-3 chain for low_relE2_w123 (relE2 features, m_a 1-3 all trained inside 115-133 GeV), so that the three
# low-mass points follow one training rule. Same steps as followup_w123.sh 3-5, with two differences:
#   - the x0.5 gate already passed in lowcand/summary_c3.txt (closure Z 0.91 / 0.47 / 0.04, data GOF min 0.24 / 0.70 / 0.32,
#     bootstrap mean/sd 0.35 / 0.09 / 0.32), so the cuts are passed directly: 0.98115 0.98837 0.99226
#   - background envelopes use the production fTest (m_a = 1 keeps Laurent, as in chainlow_relE2), no fTest_noLau_mA1 patch
# chainlow_relE2w123/ is produced beforehand by make_chainlow.py (checked line by line against chainlow_relE2).
# Steps: 1 m_a 1-3 signal re-scoring (51 condor jobs)  2 build_code_copy.sh  3 chainlow.sh under the watchdog  4 bias watch
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
L=$V/lowcand; REG=$V/logs/PENDING_RESULTS.txt
N=low_relE2_w123; TAG=relE2w123; CH=$V/chainlow_$TAG; SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
ST=$CH/logs/followup.status
stat() { echo "$(date '+%F %T') $*" > $ST; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  w123 full chain  $*" >> $REG; exit 1; }
[ -d $CH ] || stop "$CH missing (run make_chainlow.py first)"
# 1
grep -E "^mA_M[123](_|$)" $V/condor/joblist_sig.txt | awk -v o=$SBE/run3_bdt_scored_$N -v n=$N '{print $1, o, n, "high"}' > $L/list_rescore_sig_w123.txt
[ "$(wc -l < $L/list_rescore_sig_w123.txt)" -eq 51 ] || stop "signal rescore list is not 51 lines"
CL=$(cd $L && condor_submit rescore.sub LIST=list_rescore_sig_w123.txt 2>&1 | grep -oP 'submitted to cluster \K[0-9]+')
[ -n "$CL" ] || stop "signal rescore submission failed"
echo "$(date '+%F %H:%M') START w123 full chain: signal rescore cluster $CL -> $CH/logs/followup.status" >> $REG
stat "signal rescore cluster $CL"
# 2
# the copy refuses an existing target: a finished earlier build (its log ends with the du line) is kept
if [ -d $CH/flashggFinalFit ] && tail -1 $CH/logs/build_code_copy.log 2>/dev/null | grep -q "$CH/flashggFinalFit$"; then :
else bash $CH/build_code_copy.sh > $CH/logs/build_code_copy.log 2>&1 || stop "build_code_copy.sh failed ($CH/logs/build_code_copy.log)"; fi
sed -e "s#chainlow_relE2#chainlow_$TAG#g" -e "s#low_swapPt_relE2 x0.5#$N x0.5#" $V/chainlow_relE2/drive_chainlow.sh > $CH/drive_chainlow.sh
sed -e "s#chainlow_relE2#chainlow_$TAG#g" $V/chainlow_relE2/watch_bias.sh > $CH/watch_bias.sh
grep -qE "chainlow_relE2([^w]|$)" $CH/drive_chainlow.sh $CH/watch_bias.sh && stop "path substitution incomplete"
# 3 the chain only starts when every re-scored signal sample exists (drive_chainlow only waits for an empty queue)
while [ "$(condor_q $CL -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  stat "signal rescore cluster $CL: $(condor_q $CL -af ProcId 2>/dev/null | wc -l) in queue, $(condor_q $CL -hold -af ProcId 2>/dev/null | wc -l) held"; sleep 300; done
XS=$(echo $SBE/run3_bdt_scored_$N | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#')
miss=$(awk '{print $1}' $L/list_rescore_sig_w123.txt | while read s; do n=$(xrdfs eosuser.cern.ch ls $XS/$s 2>/dev/null | grep -c '\.root$'); m=$(xrdfs eosuser.cern.ch ls /eos/user/p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored_low_swapPt_relE2/$s 2>/dev/null | grep -c '\.root$'); [ "$n" -eq "$m" ] && [ "$n" -gt 0 ] || echo "$s:$n/$m"; done)
[ -z "$miss" ] || stop "signal rescore incomplete (w123/relE2 era files): $miss"
stat "chainlow running"
bash $CH/drive_chainlow.sh $CL
[ -e $CH/logs/chain4.DONE ] || stop "chainlow did not finish ($(cat $CH/logs/chain4.FAILED 2>/dev/null))"
# 4
BC=$(grep -oP 'submitted to cluster \K[0-9]+' $CH/logs/bias_submit.log | head -1)
[ -n "$BC" ] || stop "no bias cluster in bias_submit.log"
stat "bias cluster $BC"
bash $CH/watch_bias.sh $BC
stat "DONE"
