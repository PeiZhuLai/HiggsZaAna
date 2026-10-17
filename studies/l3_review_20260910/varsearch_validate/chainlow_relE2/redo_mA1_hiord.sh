#!/usr/bin/env bash
# low_swapPt_relE2, m_a = 1, option B (user 2026-10-15: "先跑B，不行就用A"): the production cap of m_a = 1 to Exp1/Pow1 is
# lifted (Laurent kept), so the envelope may add Exp3/Pow3. Without Laurent the bias moved to Exp1 (-0.267), with the
# nominal envelope Lau1 is -0.261. Redoes only m_a = 1; the current (no-Laurent) products are renamed with suffix _noLau.
# If every envelope member has |pull mean| <= 0.2 option B is kept; otherwise option A is restored automatically
# (the nominal envelope Exp1/Lau1/Pow1 from fit_results_run3/1_withLau, its limits, closure and t2w workspaces).
#   1 background envelope: sandbox fTest with the m_a = 1 order cap lifted (varsearch_validate/fTest_mA1_highorder) + makeBkgPlots
#   2 text2workspace + blind expected limits, xsec and wilson cards (cards unchanged: they reference fit_results_run3/1)
#   3 closure: pseudo-data envelope of the step-3 scan at the working point rebuilt without Laurent, same procedure as
#     chainlow.sh stage 9 (xsec card, new signal model)
#   4 bias: 10 x 1000 toys, xsec card -> EOS <MV>/bias_outputs_highstat_hiord
#   5 decision B / A (A restored as described above)
# Summary -> chainlow_relE2/redo_mA1_hiord.txt; registry DONE/STOP line.
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2
V=$(dirname $CH); S=$CH/flashggFinalFit; F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
MV=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/chainlow_relE2; FR=$S/Background/ALP_BkgModel_ReReco/fit_results_run3; DC=$S/Datacard/output_Datacard_leptons
FT=$V/fTest_mA1_highorder/bin/fTest_ALP_turnOn; L=$CH/logs/redo_mA1_hiord; mkdir -p $L
REG=$V/logs/PENDING_RESULTS.txt; ST=$L/status; OUT=$CH/redo_mA1_hiord.txt
CUT=0.98036; PSE=/eos/user/p/pelai/HZa/sandbox_varsearch_v1/closure_scan/cand_low_swapPt_relE2_drop_mA1/cut_098036
st() { echo "$(date '+%F %T') $*" > $ST; echo "$(date '+%F %T') $*" >> $L/progress.log; }
stop() { st "STOP $*"; echo "$(date '+%F %H:%M') STOP  relE2 mA1 option B  $*" >> $REG; exit 1; }
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval $(scramv1 runtime -sh)
export PYTHONPATH="$S/tools:$S/Signal/tools"
# 0 keep the previous m_a = 1 products
[ -d $FR/1_withLau ] || stop "option A backup $FR/1_withLau missing"
[ -d $FR/1_noLau ] || mv $FR/1 $FR/1_noLau || stop "cannot rename $FR/1"
for k in xsec wilson; do f=$S/Combine/output_combine_results_$k/higgsCombine1.AsymptoticLimits.mH125.38.root
  [ -e $f.withLau ] || stop "option A limit backup missing ($f.withLau)"; [ -e $f.noLau ] || cp -p $f $f.noLau; done
[ -d $CH/closure_final/mA1_noLau ] || mv $CH/closure_final/mA1 $CH/closure_final/mA1_noLau
# 1 envelope
st "1 envelope, option B"
P=$FR/1; mkdir -p $P/HZAmassInde_fTest $FR/AllFitResults_blind $FR/AllFitResults_unblind
( cd $F/Background && $FT -i $MV/root_MVAcut/data/mA_M1/ws/run3.root --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root \
    -D $P/HZAmassInde_fTest --mass_ALP 1 -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $P/fTest.dat \
    --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
for b in blind unblind; do extra=""; [ $b = unblind ] && extra="--unblind"
  ( cd $F/Background && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_$b \
      --total_OutDir $FR/AllFitResults_$b -o $P/BkgPlots_$b.root --sqrts 13p6TeV --isMultiPdf --useBinnedData \
      --massStep 2.5 --mhVal 125.0 --maVal 1 --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
      --intLumi 172.13 -c 0 --isFlashgg 0 --doBands $extra ) >> $P/makeBkgPlots.stdout.log 2>&1; done
ENV=$(grep -oE "pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+" $P/HZAmassInde_fTest/EnvelopeResults.txt | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')
[ -s $P/CMS-HGG_mva_13p6TeV_multipdf.root ] && [ -n "$ENV" ] || stop "no envelope written ($P)"
echo "$ENV" | grep -q Lau || stop "option B envelope lost Laurent: $ENV"
# 2 t2w + limits
st "2 t2w + limits"
for k in xsec wilson; do card=$DC/1_pruned_datacard_leptons.txt; [ $k = xsec ] && card=$DC/1_pruned_datacard_leptons_xsec.txt
  grep -q "$FR/1/CMS-HGG_mva_13p6TeV_multipdf.root" $card || stop "$k card does not reference $FR/1"
  mkdir -p $S/Combine/work_$k/1
  ( cd $S/Combine/work_$k/1 && text2workspace.py $card -o $S/Combine/root_t2w_$k/1_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > $L/t2w_$k.log 2>&1 \
    && combine -M AsymptoticLimits $S/Combine/root_t2w_$k/1_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
         --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n 1 > $L/limit_$k.log 2>&1 \
    && mv higgsCombine1.AsymptoticLimits.mH125.38.root $S/Combine/output_combine_results_$k/ ) || stop "t2w/limit $k failed ($L)"
done
# 3 closure
st "3 closure"
PSD=$CH/closure_final/mA1_psfit_hiord; O=$CH/closure_final/mA1; mkdir -p $PSD/ft $O
xrdcp -f -s root://eosuser.cern.ch/$PSE/ps/mA_M1/ws/run3.root $PSD/ps_ws_run3.root || stop "cannot fetch the pseudo-data workspace"
( cd $F/Background && $FT -i $PSD/ps_ws_run3.root --saveMultiPdf $PSD/CMS-HGG_mva_13p6TeV_multipdf.root -D $PSD/ft --mass_ALP 1 -c 1 \
    --isFlashgg 0 --isData 0 -f data, --datfile $PSD/fTest.dat --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $PSD/ftest.stdout.log 2>&1
PSENV=$(grep -oE "pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+" $PSD/ft/EnvelopeResults.txt | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')
[ -s $PSD/CMS-HGG_mva_13p6TeV_multipdf.root ] || stop "no pseudo-data envelope"
sed "s#$FR/1/#$PSD/#g" $DC/1_pruned_datacard_leptons_xsec.txt > $O/datacard.txt
grep -q "$PSD/CMS-HGG_mva_13p6TeV_multipdf.root" $O/datacard.txt || stop "closure: background path not swapped"
( cd $O && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
  combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n PS1 > sig.log 2>&1
  combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n PS1 > rhat.log 2>&1
  rm -f ws_t2w.root )
CLO="closure Z=$(grep -E 'Significance:' $O/sig.log | sed -E 's/.*Significance: //')  r_hat=$(grep -E '^\s+r :' $O/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')"
# 4 bias
st "4 bias submit"
HS_OUT_BASE=$MV/bias_outputs_highstat_hiord ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec NCHUNKS=10 BASE_DIR=$S \
  bash $S/shellScripts/bias/Condor/subjob_bias_highstat.sh 1 > $L/bias_submit.log 2>&1
BC=$(grep -oP 'submitted to cluster \K[0-9]+' $L/bias_submit.log | head -1); [ -n "$BC" ] || stop "bias submission failed ($L/bias_submit.log)"
while [ "$(condor_q $BC -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  s=$(condor_q $BC -af JobStatus 2>/dev/null); h=$(echo "$s" | grep -c '^5$'); r=$(echo "$s" | grep -c '^2$'); t=$(echo "$s" | grep -c .)
  st "4 bias cluster $BC: $t in queue, $r running, $h held"
  [ "$t" -gt 0 ] && [ "$h" -eq "$t" ] && stop "all $t remaining bias jobs held"
  sleep 600; done
HS_OUT_BASE=$MV/bias_outputs_highstat_hiord ROOT_DATACARD_PATH=$S/Combine/root_t2w_xsec BASE_DIR=$S \
  bash $S/shellScripts/bias/Condor/merge_bias_highstat.sh 1 > $L/bias_merge.log 2>&1
python3 - $S/Combine $MV "$ENV" "$PSENV" "$CLO" > $OUT 2>&1 <<'PY' || stop "summary failed ($OUT)"
import json, sys, ROOT
ROOT.gROOT.SetBatch(True)
C, MV, env, psenv, clo = sys.argv[1:6]
def q(p):
    f = ROOT.TFile.Open(p); o = {round(e.quantileExpected, 3): e.limit for e in f.Get("limit")}; f.Close(); return o
print("low_swapPt_relE2, m_a = 1 at cut 0.98036, option B: Exp/Pow order cap lifted, Laurent kept")
print("data envelope:", env)
for k in ("xsec", "wilson"):
    n = "%s/output_combine_results_%s/higgsCombine1.AsymptoticLimits.mH125.38.root" % (C, k)
    a, b = q(n), q(n + ".withLau")   # b = option A
    print("expected limit %-6s: %.4f [%.4f, %.4f]  (with Laurent %.4f, ratio %.3f)" % (k, a[0.5], a[0.16], a[0.84], b[0.5], a[0.5] / b[0.5]))
print("closure (pseudo-data envelope, option B: %s): %s" % (psenv, clo))
def bias(base):
    try: return json.load(open("%s/mA_1/merged/BiasJson/1_gaussfit.json" % base))["fit_results"]
    except OSError: return None
for tag, base in (("option B", MV + "/bias_outputs_highstat_hiord"), ("A: Exp1/Lau1/Pow1", MV + "/bias_outputs_highstat"), ("without Laurent", MV + "/bias_outputs_highstat_noLau")):
    d = bias(base)
    print("bias pull mean, %-20s:" % tag, " ".join("%s:%+.3f%s" % (k, v["mean"], "!" if abs(v["mean"]) > 0.2 else "") for k, v in sorted(d.items())) if d else "MISSING")
PY
grep -q "MISSING" $OUT && stop "bias result missing ($OUT)"
# 5 decision
worst=$(python3 -c "
import json; d = json.load(open('$MV/bias_outputs_highstat_hiord/mA_1/merged/BiasJson/1_gaussfit.json'))['fit_results']
print('%.3f' % max(abs(v['mean']) for v in d.values()))")
if python3 -c "import sys; sys.exit(0 if $worst <= 0.2 else 1)"; then
  echo "DECISION: option B kept (max |pull mean| $worst <= 0.2)" >> $OUT; dec="B kept"
else
  st "5 restore option A"
  mv $FR/1 $FR/1_hiord && cp -rp $FR/1_withLau $FR/1 || stop "cannot restore the option A envelope"
  mv $CH/closure_final/mA1 $CH/closure_final/mA1_hiord && cp -rp $CH/closure_final/mA1_withLau $CH/closure_final/mA1
  for k in xsec wilson; do card=$DC/1_pruned_datacard_leptons.txt; [ $k = xsec ] && card=$DC/1_pruned_datacard_leptons_xsec.txt
    f=$S/Combine/output_combine_results_$k/higgsCombine1.AsymptoticLimits.mH125.38.root; cp -p $f $f.hiord
    ( cd $S/Combine/work_$k/1 && text2workspace.py $card -o $S/Combine/root_t2w_$k/1_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > $L/t2w_A_$k.log 2>&1 \
      && combine -M AsymptoticLimits $S/Combine/root_t2w_$k/1_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
           --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n 1 > $L/limit_A_$k.log 2>&1 \
      && mv higgsCombine1.AsymptoticLimits.mH125.38.root $f ) || stop "restoring option A: t2w/limit $k failed"
  done
  python3 - $S/Combine >> $OUT <<'PYA' || stop "option A not reproduced"
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def q(p):
    f = ROOT.TFile.Open(p); o = {round(e.quantileExpected, 3): e.limit for e in f.Get("limit")}; f.Close(); return o
for k in ("xsec", "wilson"):
    n = "%s/output_combine_results_%s/higgsCombine1.AsymptoticLimits.mH125.38.root" % (sys.argv[1], k)
    a, b = q(n), q(n + ".withLau")
    print("restored option A %s limit %.4f (previous option A %.4f)" % (k, a[0.5], b[0.5]))
    assert abs(a[0.5] / b[0.5] - 1) < 1e-6
PYA
  echo "DECISION: option B not kept (max |pull mean| $worst > 0.2) -> option A restored (envelope, t2w, limits, closure; option B products kept as *_hiord)" >> $OUT; dec="A restored"
fi
st "DONE ($dec)"
echo "$(date '+%F %H:%M') DONE  relE2 mA1 option B -> $dec  -> $OUT" >> $REG
