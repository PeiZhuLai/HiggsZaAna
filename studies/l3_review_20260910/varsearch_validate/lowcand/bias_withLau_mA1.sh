#!/usr/bin/env bash
# Cross-check for the m_a = 1 envelope of low_relE2_w123 (side-agent note accepted by the user 2026-10-15): the nominal
# chain drops Laurent at m_a = 1 (Lau1 -0.261 measured on the previous low_swapPt_relE2 model, window training at
# m_a = 1 only). This builds the m_a = 1 envelope of the final model with the production fTest (Laurent kept), a copy of
# the xsec card pointing at it, and runs the same high-statistics bias (10 x 1000 toys) so that the decision rests on
# the final model's own Lau1 number.
# Waits for chainlow_relE2w123/logs/chain4.DONE (stops if the follow-up driver stopped). Outputs:
#   flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3_withLau/1/   envelope with Laurent
#   flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons_xsec_withLau.txt, Combine/root_t2w_xsec_withLau/
#   EOS <MV>/bias_outputs_highstat_withLau, summary chainlow_relE2w123/bias_withLau_mA1.txt; registry line.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
CH=$V/chainlow_relE2w123; S=$CH/flashggFinalFit; MV=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2w123
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
REG=$V/logs/PENDING_RESULTS.txt; ST=$V/lowcand/logs/bias_withLau.status
st() { echo "$(date '+%F %T') $*" > $ST; }
stop() { st "STOP $*"; echo "$(date '+%F %H:%M') STOP  w123 mA1 bias with Laurent  $*" >> $REG; exit 1; }
while [ ! -e $CH/logs/chain4.DONE ]; do
  grep -q STOP $V/lowcand/logs/followup_w123.status 2>/dev/null && stop "follow-up driver stopped"
  st "waiting for chainlow_relE2w123"; sleep 600; done
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval $(scramv1 runtime -sh)
FR=$S/Background/ALP_BkgModel_ReReco/fit_results_run3; A=$S/Background/ALP_BkgModel_ReReco/fit_results_run3_withLau/1
mkdir -p $A/HZAmassInde_fTest
st "fTest with Laurent"
( cd $F/Background && ./bin/fTest_ALP_turnOn -i $MV/root_MVAcut/data/mA_M1/ws/run3.root --saveMultiPdf $A/CMS-HGG_mva_13p6TeV_multipdf.root \
    -D $A/HZAmassInde_fTest --mass_ALP 1 -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $A/fTest.dat \
    --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $A/ftest.stdout.log 2>&1
grep -q "Lau" $A/HZAmassInde_fTest/EnvelopeResults.txt || stop "envelope with Laurent has no Lau member ($A)"
D=$S/Datacard/output_Datacard_leptons; X=$D/1_pruned_datacard_leptons_xsec.txt; XL=$D/1_pruned_datacard_leptons_xsec_withLau.txt
[ "$(grep -c "$FR/1/CMS-HGG_mva_13p6TeV_multipdf.root" $X)" -eq 2 ] || stop "xsec card does not reference the m_a 1 envelope twice"
sed "s#$FR/1/CMS-HGG_mva_13p6TeV_multipdf.root#$A/CMS-HGG_mva_13p6TeV_multipdf.root#g" $X > $XL
mkdir -p $S/Combine/root_t2w_xsec_withLau $S/Combine/work_withLau
( cd $S/Combine/work_withLau && text2workspace.py $XL -o $S/Combine/root_t2w_xsec_withLau/1_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 ) > $CH/logs/t2w_withLau.log 2>&1 \
  || stop "text2workspace failed ($CH/logs/t2w_withLau.log)"
st "bias submit"
HS_OUT_BASE=$MV/bias_outputs_highstat_withLau ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec_withLau NCHUNKS=10 \
  bash $S/shellScripts/bias/Condor/subjob_bias_highstat.sh 1 > $CH/logs/bias_withLau_submit.log 2>&1
CL=$(grep -oP 'submitted to cluster \K[0-9]+' $CH/logs/bias_withLau_submit.log | head -1); [ -n "$CL" ] || stop "bias submission failed"
while [ "$(condor_q $CL -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  s=$(condor_q $CL -af JobStatus 2>/dev/null); h=$(echo "$s" | grep -c '^5$'); t=$(echo "$s" | grep -c .)
  st "bias cluster $CL in queue $t held $h"; [ "$t" -gt 0 ] && [ "$h" -eq "$t" ] && stop "all $t remaining bias jobs held"
  sleep 600; done
HS_OUT_BASE=$MV/bias_outputs_highstat_withLau ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec_withLau BASE_DIR=$S \
  bash $S/shellScripts/bias/Condor/merge_bias_highstat.sh 1 > $CH/logs/bias_withLau_merge.log 2>&1
python3 - $MV/bias_outputs_highstat_withLau $A/HZAmassInde_fTest/EnvelopeResults.txt > $CH/bias_withLau_mA1.txt 2>&1 <<'PY' || stop "no merged bias result"
import json, re, sys
d = json.load(open("%s/mA_1/merged/BiasJson/1_gaussfit.json" % sys.argv[1]))["fit_results"]
print("m_a = 1, low_relE2_w123, envelope WITH Laurent (production fTest), xsec card, 10 x 1000 toys")
print("envelope:", " ".join("%s:%s" % e for e in re.findall(r"pdf : ([A-Za-z0-9]+) , gof : ([0-9.e-]+)", open(sys.argv[2]).read())))
print("pull mean per generating function:", " ".join("%s:%+.3f%s" % (k, v["mean"], "!" if abs(v["mean"]) > 0.2 else "") for k, v in sorted(d.items())))
PY
st "DONE"
echo "$(date '+%F %H:%M') DONE  w123 mA1 bias with Laurent  -> $CH/bias_withLau_mA1.txt" >> $REG
