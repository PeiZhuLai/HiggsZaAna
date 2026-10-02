#!/usr/bin/env bash
# 2026-10-01 20:0x: mA3 datacard did not carry the S1 scheme-2 lnN. chain_dyveto_wp3_env2425.sh re-ran
# apply_bdt_sig.py for mA3 without HZA_MVA_REWEIGHT_LNN_JSON (-> old per-event variation), and the lnN
# JSON had been derived at the old mA3 cut 0.992. Here:
#   1. S1 evaluation for mA3 at 0.988 -> replace the mA3 entry of s1_lnn_scheme2_fsrfix_dyveto.json
#   2. apply_bdt_sig.py mA3 (5 eras) WITH the JSON; check the constants in every signal tree
#   3. Tree2WS mA3 (signal fits are unaffected: the lnN only changes the weight columns)
#   4. all datacards (check the mA3 row) + t2w + limits
#   5. condor: bias + impact mA3, impact mA29 (its impacts json had only 51 parameters); merge mA3 bias
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FZ=$C/root_t2w_dyveto_20260930
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
S1D=$D/s1_altsideband_20260927
LNN=$S1D/results/s1_lnn_scheme2_fsrfix_dyveto.json
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
TAG=before_mA3_s1fix_$(date +%Y%m%d_%H%M)
HB=$L/dyveto_mA3_s1fix_heartbeat.txt
S1=$D/wp3_single; mkdir -p $S1
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }
die(){ hb "STOPPING: $*"; exit 1; }
fresh(){ [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

hb "[1] S1 evaluation mA3 at the current working point"
grep -q '"mA": 3, "MVAcut": 0.988\|"MVAcut": 0.988' /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json || \
  $ALP_PY -c "import json,sys; r={e['mA']:e['MVAcut'] for e in json.load(open('/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json'))['results']}; sys.exit(0 if abs(r[3]-0.988)<1e-6 else 1)" || die "MVAcut JSON mA3 is not 0.988"
env -i HOME=$HOME PATH=/usr/bin:/bin S1_RW_TAG=fsrfix_dyveto $ALP_PY $S1D/eval_alt_reweight.py --masses 3 --out $S1D/results/s1_altsideband_yields_dyveto_mA3_0988.csv > $L/s1fix_eval_mA3.log 2>&1 || die "S1 eval mA3"
env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY $S1D/make_s1_lnn_json.py $S1D/results/s1_altsideband_yields_dyveto_mA3_0988.csv $S1D/results/s1_lnn_mA3_0988.json > $L/s1fix_lnn_mA3.log 2>&1 || die "lnN mA3"
cat $L/s1fix_lnn_mA3.log | tee -a $HB
cp -n $LNN $LNN.$TAG
$ALP_PY - $LNN $S1D/results/s1_lnn_mA3_0988.json <<'PY' || die "merge mA3 lnN"
import json, sys
main, m3 = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
old = main["values"]["3"]; main["values"]["3"] = m3["values"]["3"]
main["note_mA3"] = "mA3 re-evaluated at the 0.988 working point (2026-10-01); was %s at 0.992" % old
json.dump(main, open(sys.argv[1], "w"), indent=1)
print("mA3 lnN %s -> %s" % (old, main["values"]["3"]))
PY
export HZA_MVA_REWEIGHT_LNN_JSON=$LNN

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"

hb "[2] apply_bdt_sig mA3 with the lnN JSON"
for y in $YEARS; do cp -n $EOSMVA/sig/mA_M3/output_$y.root $EOSMVA/sig/mA_M3/output_$y.root.$TAG 2>/dev/null; done
T=$(date +%s)
for y in $YEARS; do python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M3 --years $y > $L/s1fix_mva_sig_$y.log 2>&1 & done
wait
for y in $YEARS; do fresh $EOSMVA/sig/mA_M3/output_$y.root $T || die "signal MVA cut mA3 $y stale"; done
python3 - $EOSMVA $LNN <<'PY' || die "lnN not in mA3 signal trees"
import json, sys, uproot, numpy as np
B, J = sys.argv[1], json.load(open(sys.argv[2]))["values"]["3"]
for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
    f = uproot.open("%s/sig/mA_M3/output_%s.root" % (B, y))
    for lep in ("ele", "mu"):
        a = f["DiphotonTree/ggh_125_Za_%s_13p6TeV_cat0" % lep].arrays(["weight_mva_reweight_Up", "weight_mva_reweight_Down"], library="np")
        ok = np.allclose(a["weight_mva_reweight_Up"], J["up"], atol=1e-5) and np.allclose(a["weight_mva_reweight_Down"], J["down"], atol=1e-5)
        print("  mA3 %s %s up/down %s" % (y, lep, "OK" if ok else "MISMATCH"))
        if not ok: sys.exit(1)
PY

hb "[3] Tree2WS mA3"
sed -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(3)/' \
    -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(3)/' \
    $F/Trees2WS/run_tree2ws.sh > $S1/run_tree2ws_m3.sh
grep -q '^mAs_sig=(3)$' $S1/run_tree2ws_m3.sh || die "tree2ws single-mass copy failed"
T=$(date +%s)
( cd $F/Trees2WS && bash $S1/run_tree2ws_m3.sh > $L/s1fix_tree2ws.log 2>&1 )
n=$(find $EOSMVA/sig/mA_M3/ws_Tree2WS -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 10 ] || die "only $n fresh signal ws for mA3"

hb "[4] all datacards + t2w + limits"
V=$C/output_combine_results_variants; mkdir -p $V
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
cp -rp $F/Datacard/output_Datacard_leptons $F/Datacard/output_Datacard_leptons_$TAG
T=$(date +%s)
( cd $F/Datacard && sh 1_runjob_gen_datacard_makeYields.sh > $L/s1fix_dc_1.log 2>&1 && sh 2_runjob_gen_datacard_makeDatacard.sh > $L/s1fix_dc_2.log 2>&1 && sh 3_rysn_datacard.sh > $L/s1fix_dc_3.log 2>&1 )
for m in $(seq 1 30); do fresh $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt $T || die "datacard mA$m stale"; done
row=$(grep "^CMS_hza_mva_reweight" $F/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt)
want=$(python3 -c "import json; v=json.load(open('$LNN'))['values']['3']; print('%.3f/%.3f' % (v['down'], v['up']))")
echo "$row" | grep -q "$want" || die "mA3 datacard row ($row) does not carry $want"
hb "  mA3 datacard: $(echo $row | awk '{print $1,$2,$3,$4}')"
T=$(date +%s)
( cd $C && sh 2_text2ws.sh > $L/s1fix_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/s1fix_limits.log 2>&1 )
for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit mA$m stale"; done
python3 - $C/output_combine_results $V/output_combine_results_$TAG <<'PY' | tee -a $HB
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); v = None
    for e in f.Get("limit"):
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    return v
for m in range(1, 31):
    n, o = med(sys.argv[1], m), med(sys.argv[2], m)
    if abs(n / o - 1) > 1e-4: print("  limit moved mA%d %.5f -> %.5f (x%.4f)" % (m, o, n, n / o))
print("  limit comparison done")
PY

hb "[5] condor: bias + impact mA3, impact mA29"
cp -p $FZ/3_Datacard_leptons.root $FZ/3_Datacard_leptons.root.$TAG 2>/dev/null
cp -p $C/root_t2w/3_Datacard_leptons.root $FZ/
cp -p $FZ/29_Datacard_leptons.root $FZ/29_Datacard_leptons.root.$TAG 2>/dev/null
cp -p $C/root_t2w/29_Datacard_leptons.root $FZ/
mkdir -p $BN/bias_outputs_highstat_variants
[ -d $BN/bias_outputs_highstat/mA_3 ] && mv $BN/bias_outputs_highstat/mA_3 $BN/bias_outputs_highstat_variants/mA_3_$TAG
cp -p $C/output_impacts/29_impacts.json $C/output_impacts/29_impacts.json.$TAG 2>/dev/null
T0=$(date +%s)
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh 3 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh 3 29 2>&1)
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "submission failed"
hb "  clusters $CL"
while true; do
  q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  while read -r id mem; do
    [ -z "${id:-}" ] && continue
    if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
    hb "  held $id mem $mem -> $new"; condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null
  done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
  hb "  idle=$idle running=$run held=$held"
  [ $((idle + run + held)) -eq 0 ] && break
  [ $((idle + run)) -eq 0 ] && [ "$held" -gt 0 ] && { hb "  only held left"; break; }
  sleep 1800
done
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh 3 > $L/s1fix_bias_merge.log 2>&1
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/s1fix_bias_collect.log 2>&1
python3 - $BN/bias_outputs_highstat $C/output_impacts $T0 <<'PY' | tee -a $HB
import json, os, sys
B, I, T0 = sys.argv[1], sys.argv[2], float(sys.argv[3]); bad = 0
p = "%s/mA_3/merged/BiasJson/3_gaussfit.json" % B
if os.path.exists(p):
    fr = json.load(open(p))["fit_results"]
    items = sorted(((k, v["mean"]) for k, v in fr.items()), key=lambda x: -abs(x[1]))
    print("  mA3 bias " + " ".join("%s:%+.3f" % x for x in items)); bad += any(abs(v) > 0.2 for _, v in items)
else:
    print("  mA3 bias MISSING"); bad += 1
for m in (3, 29):
    ip = "%s/%d_impacts.json" % (I, m)
    if os.path.exists(ip) and os.path.getmtime(ip) >= T0:
        d = json.load(open(ip)); r = [q for q in d["POIs"] if q["name"] == "r"][0]["fit"][1]
        print("  mA%d impact r=%.3f, %d parameters" % (m, r, len(d["params"])))
    else:
        print("  mA%d impact MISSING/stale" % m); bad += 1
print("VERDICT: %s" % ("PASSED" if bad == 0 else "CHECK"))
PY
echo "$(date '+%F %H:%M') DONE mA3 S1 fix + mA29 impact (heartbeat dyveto_mA3_s1fix_heartbeat.txt)" >> $L/PENDING_RESULTS.txt
