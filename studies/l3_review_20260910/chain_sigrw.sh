#!/usr/bin/env bash
# Signal reweighting (user decision 2026-10-02): the signal is reweighted with the nominal sideband
# (data/MC) reweight, like the background simulation, normalized per (m_a, era, channel) to keep the
# preselected yield; CMS_hza_mva_reweight = envelope of {z_m&hsb, z_m_narrow&hsb, no reweight} relative
# to the nominal reweight (option B), s1_lnn_scheme2_sigrw_dyveto.json.
# Code: apply_bdt_sig.py (signal reweight + reference check), signal_eff_sumw.py (same weights for the
# interpolation curves), make_s1_lnn_json.py --reference nominal. Validated on mA5 2024: efficiency equal
# to eval_alt_reweight eff_nominal to 1e-16; signal_eff_sumw equal to apply_bdt_sig (muon exact).
# Stages (START_STAGE=N skips < N):
#   1 snapshot / park previous products        5 signal model (140 fits)
#   2 interpolation efficiency JSON             6 datacards, t2w, limits (compare), limit plots
#   3 MVA cut, signal only (70 files) + checks  7 condor bias (30x10) + impacts (30); closure locally
#   4 Tree2WS, signal only                      8 summary
# Data MVA cut / workspaces and background envelopes are unchanged (data and background untouched).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
PC=$F/pseudodata_closure
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
LNN=$D/s1_altsideband_20260927/results/s1_lnn_scheme2_sigrw_dyveto.json
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
TAG=preSigRW_20261002
SNAP=$D/snapshot_$TAG
FZ=$C/root_t2w_sigrw_20261002
HB=$L/sigrw_heartbeat.txt
MASS14="1 2 3 4 5 6 7 8 9 10 15 20 25 30"
MASS30=$(seq 1 30 | tr '\n' ' ')
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
EFFJ_ELE=$PLOT/output/sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json
EFFJ_MU=$PLOT/output/sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json
START=${START_STAGE:-1}
run() { [ "$START" -le "$1" ]; }
hb() { echo "$(date '+%F %T') $*" | tee -a "$HB"; }
die() { hb "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

if run 1; then
hb "[1] snapshot ($SNAP) and park previous signal products"
mkdir -p $SNAP
cp -rp $C/output_combine_results $SNAP/
cp -rp $F/Datacard/output_Datacard_leptons $SNAP/
for ch in ele mu; do mkdir -p $SNAP/sigfit_$ch && cp -p $F/Signal/outdir_$ch/signalFit/output/*.root $SNAP/sigfit_$ch/ 2>/dev/null; done
cp -p $EFFJ_ELE $EFFJ_ELE.$TAG; cp -p $EFFJ_MU $EFFJ_MU.$TAG
[ -d $EOSMVA/sig_$TAG ] || mv $EOSMVA/sig $EOSMVA/sig_$TAG
mkdir -p $EOSMVA/sig
mkdir -p $PC/$TAG && for d in closure_results.json output_limits output_significance output_fitdiag plots; do [ -e $PC/$d ] && cp -rp $PC/$d $PC/$TAG/; done
hb "    done"
fi

if run 2; then
hb "[2] interpolation efficiency JSON (signal_eff_sumw.py, reweighted signal)"
T=$(date +%s)
( cd $PLOT && env -u PYTHONHOME PYTHONPATH=$PLOT/lib $ALP_PY scripts/signal_eff_sumw.py > $L/sigrw_signal_eff_sumw.log 2>&1 )
fresh $EFFJ_ELE $T && fresh $EFFJ_MU $T || die "efficiency JSON not rewritten (sigrw_signal_eff_sumw.log)"
env -u PYTHONPATH -u PYTHONHOME $ALP_PY - "$EFFJ_ELE" "$EFFJ_MU" "$EFFJ_ELE.$TAG" <<'PY' | tee -a $HB || die "efficiency JSON check failed"
import json, sys
new_e, new_m, old_e = (json.load(open(p)) for p in sys.argv[1:4])
for d, n in ((new_e, "ele"), (new_m, "mu")):
    v = d["values"]; nz = sum(1 for y in v for m in range(1, 31) if float(v[y].get(str(m), 0)) > 0)
    if len(v) != 5 or nz != 150: print("  %s: %d years %d/150 non-zero" % (n, len(v), nz)); sys.exit(1)
r = [float(new_e["values"][y][str(m)]) / float(old_e["values"][y][str(m)]) for y in new_e["values"] for m in range(1, 31)]
print("    ele efficiency new/old: %.3f - %.3f" % (min(r), max(r)))
PY
fi

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"
set -u

if run 3; then
hb "[3] MVA cut, signal (14 mA x 5 eras) with the signal reweight"
# Resume after a partial failure: STAGE3_T0=<epoch of the original stage-3 start> SIG_APPLY_ONLY="mA_M10:2022preEE ..."
# reruns only those jobs; the freshness check below still uses the original start time.
T=${STAGE3_T0:-$(date +%s)}
SCR=$F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py
ALOG=$F/MVAcut/run3_ReReco_Sys/logs/apply_bdt_sig
if [ -n "${SIG_APPLY_ONLY:-}" ]; then
  failed="$SIG_APPLY_ONLY"
else
  bash $F/shellScripts/mva/run_apply_bdt_sig_6jobs.sh > $L/sigrw_apply_sig.log 2>&1
  failed=$(grep -oP '^\[ERROR\] Failed \K(mA_M[0-9]+)_(\S+?)(?=;)' $L/sigrw_apply_sig.log | sed -E 's/^(mA_M[0-9]+)_/\1:/' | tr '\n' ' ')
fi
# transient EOS I/O errors: rerun each failed job on its own, up to 3 times
for job in $failed; do
  ok=0
  for att in 1 2 3; do
    hb "  rerun ${job%%:*} ${job##*:} (attempt $att)"
    if python3 $SCR --samples ${job%%:*} --years ${job##*:} > $ALOG/${job%%:*}_${job##*:}.log 2>&1; then ok=1; break; fi
    sleep 60
  done
  [ $ok -eq 1 ] || die "apply_bdt_sig ${job} failed 3 times ($ALOG/${job%%:*}_${job##*:}.log)"
done
bad=0
for m in $MASS14; do for y in $YEARS; do fresh $EOSMVA/sig/mA_M$m/output_$y.root $T || { echo "  sig mA_M$m $y stale/missing"; bad=1; }; done; done
[ $bad -eq 0 ] || die "signal MVA cut outputs incomplete (sigrw_apply_sig.log, MVAcut/run3_ReReco_Sys/logs/apply_bdt_sig)"
python3 $D/audit_mvacut_trees.py > $L/sigrw_mva_audit.log 2>&1 || { tail -15 $L/sigrw_mva_audit.log; die "MVA cut trees fail a full read"; }
python3 - "$EOSMVA" "$LNN" <<'PY' | tee -a $HB || die "signal trees lack the reweight or the new lnN"
import json, sys, uproot, numpy as np
B, J = sys.argv[1], json.load(open(sys.argv[2]))["values"]
bad = 0; rmin, rmax = 9, 0
for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30):
    for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        f = uproot.open("%s/sig/mA_M%d/output_%s.root" % (B, m, y))
        trees = sorted(set(k.split(";")[0] for k in f.keys() if "ggh_125_Za" in k))
        if len(trees) != 34: bad += 1; print("  mA%d %s: %d trees" % (m, y, len(trees)))
        for k in trees:
            t = f[k]
            if "weight_sideband_rwgt" not in t.keys(): bad += 1; print("  no rwgt", m, y, k); continue
        for lep in ("ele", "mu"):
            a = f["DiphotonTree/ggh_125_Za_%s_13p6TeV_cat0" % lep].arrays(
                ["weight_mva_reweight_Up", "weight_mva_reweight_Down", "weight_sideband_rwgt"], library="np")
            if len(a["weight_sideband_rwgt"]):
                rmin = min(rmin, a["weight_sideband_rwgt"].mean()); rmax = max(rmax, a["weight_sideband_rwgt"].mean())
            if len(a["weight_mva_reweight_Up"]) and not (np.allclose(a["weight_mva_reweight_Up"], J[str(m)]["up"], atol=1e-5)
                                                       and np.allclose(a["weight_mva_reweight_Down"], J[str(m)]["down"], atol=1e-5)):
                bad += 1; print("  lnN mismatch mA%d %s %s" % (m, y, lep))
print("    signal trees: %d problems; mean selected-event reweight %.3f - %.3f" % (bad, rmin, rmax))
sys.exit(1 if bad else 0)
PY
fi

if run 4; then
hb "[4] Tree2WS, signal only"
T=$(date +%s)
mkdir -p $D/wp3_single
sed -e 's/^mAs_data=(.*)$/mAs_data=()/' $F/Trees2WS/run_tree2ws.sh > $D/wp3_single/run_tree2ws_sigonly.sh
grep -q '^mAs_data=()$' $D/wp3_single/run_tree2ws_sigonly.sh || die "signal-only Tree2WS copy failed"
( cd $F/Trees2WS && bash $D/wp3_single/run_tree2ws_sigonly.sh > $L/sigrw_tree2ws.log 2>&1 )
bad=0
for m in $MASS14; do
  n=$(find $EOSMVA/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] || { echo "  sig ws mA_M$m: only $n fresh"; bad=1; }
done
[ $bad -eq 0 ] || die "signal workspaces incomplete (sigrw_tree2ws.log)"
hb "    signal workspaces fresh for 14 mass points"
fi

if run 5; then
hb "[5] signal model (fTest -> photon syst -> signalFit -> plotter -> effSigma)"
T=$(date +%s)
for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter 5_runjob_sig_plotEffSigma; do
  hb "    $s"
  ( cd $F/shellScripts && bash $F/shellScripts/sig_sys/$s.sh > $L/sigrw_sig_$s.log 2>&1 )
done
bad=0
for m in $MASS14; do for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || { echo "  sigfit mA$m $y $ch stale"; bad=1; }
done; done; done
[ $bad -eq 0 ] || die "signal fits incomplete (sigrw_sig_*.log)"
hb "    signal fits 140/140 fresh"
fi

if run 6; then
hb "[6] datacards, t2w, limits"
T=$(date +%s)
V=$C/output_combine_results_variants; mkdir -p $V
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
( cd $F/Datacard && sh 1_runjob_gen_datacard_makeYields.sh > $L/sigrw_dc_1.log 2>&1 && sh 2_runjob_gen_datacard_makeDatacard.sh > $L/sigrw_dc_2.log 2>&1 && sh 3_rysn_datacard.sh > $L/sigrw_dc_3.log 2>&1 )
bad=0; for m in $MASS30; do fresh $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt $T || { echo "  datacard mA$m stale"; bad=1; }; done
for m in $MASS30; do
  n=$(grep -c 'rateParam' $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt)
  case " $MASS14 " in *" $m "*) want=0 ;; *) want=10 ;; esac
  [ "$n" -eq "$want" ] || { echo "  datacard mA$m: $n rateParam lines (want $want)"; bad=1; }
done
[ $bad -eq 0 ] || die "datacards incomplete (sigrw_dc_*.log)"
want=$($ALP_PY -c "import json; v=json.load(open('$LNN'))['values']['1']; print('%.3f/%.3f' % (v['down'], v['up']))")
grep "^CMS_hza_mva_reweight" $F/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt | grep -q "$want" || die "mA1 datacard does not carry the new lnN $want"
T=$(date +%s)
( cd $C && sh 2_text2ws.sh > $L/sigrw_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/sigrw_limits.log 2>&1 )
bad=0; for m in $MASS30; do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || { echo "  limit mA$m stale"; bad=1; }; done
[ $bad -eq 0 ] || die "limits incomplete (sigrw_limits.log)"
python3 - "$C/output_combine_results" "$SNAP/output_combine_results" <<'PY' | tee $L/sigrw_limits_compare.txt | tee -a $HB
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); v = None
    for e in f.Get("limit"):
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
print("    mA  median_old  median_new  new/old")
for m in range(1, 31):
    o, n = med(sys.argv[2], m), med(sys.argv[1], m)
    print("    %2d  %10.5f  %10.5f  %7.3f" % (m, o, n, n / o))
PY
( cd $F/Plots && sh 2_runLimitsPlot.sh > $L/sigrw_limit_plots.log 2>&1 )
fi

if run 7; then
hb "[7] condor bias (30 x 10 chunks) + impacts (30); pseudo-data closure locally"
T0=$(date +%s)
mkdir -p $FZ && cp -p $C/root_t2w/*_Datacard_leptons.root $FZ/
for m in $MASS30; do fresh $FZ/${m}_Datacard_leptons.root $(date -d '-3 hours' +%s) || die "frozen workspace mA$m is old"; done
[ -d $BN/bias_outputs_highstat ] && [ ! -d $BN/bias_outputs_highstat_$TAG ] && mv $BN/bias_outputs_highstat $BN/bias_outputs_highstat_$TAG
[ -d $BN/plots_bias ] && [ ! -d $BN/plots_bias_$TAG ] && mv $BN/plots_bias $BN/plots_bias_$TAG
[ -d $C/output_impacts ] && [ ! -d $C/output_impacts_$TAG ] && cp -rp $C/output_impacts $C/output_impacts_$TAG
bout=$(ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $MASS30 2>&1)
iout=$(bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $MASS30 2>&1)
echo "$bout" > $L/sigrw_bias_submit.log; echo "$iout" > $L/sigrw_impact_submit.log
CL=$(grep -ohP 'submitted to cluster \K[0-9]+' <<< "$bout
$iout" | sort -u | tr '\n' ' ')
[ -n "${CL// /}" ] || die "nothing submitted (sigrw_*_submit.log)"
hb "    clusters $CL"
# closure (background-only pseudo-data is unchanged; only the fits with the new signal model rerun)
( cd $PC/scripts && START=3 bash chain_closure30_dyveto.sh > $L/sigrw_closure.log 2>&1 ) &
cpid=$!
while true; do
  q=$(condor_q $CL -af JobStatus 2>/dev/null) || { sleep 300; continue; }
  idle=$(grep -c '^1$' <<<"$q"); run_=$(grep -c '^2$' <<<"$q"); held=$(grep -c '^5$' <<<"$q")
  while read -r id mem; do
    [ -z "${id:-}" ] && continue
    if [ "$mem" -lt 8000 ]; then new=8000; elif [ "$mem" -lt 16000 ]; then new=16000; elif [ "$mem" -lt 20000 ]; then new=20000; else continue; fi
    hb "    held $id mem $mem -> $new"; condor_qedit $id RequestMemory $new >/dev/null; condor_release $id >/dev/null
  done < <(condor_q $CL -constraint 'JobStatus==5 && regexp("emory", HoldReason)' -af:j RequestMemory 2>/dev/null)
  imp=0; for m in $MASS30; do fresh $C/output_impacts/${m}_impacts.json $T0 && imp=$((imp+1)); done
  hb "    idle=$idle running=$run_ held=$held | impacts $imp/30 | closure $(kill -0 $cpid 2>/dev/null && echo running || echo done)"
  [ $((idle + run_ + held)) -eq 0 ] && break
  [ $((idle + run_)) -eq 0 ] && [ "$held" -gt 0 ] && { hb "    only held jobs left"; break; }
  sleep 1800
done
wait $cpid; hb "    closure exit $?"
ROOT_DATACARD_PATH=$FZ bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $MASS30 > $L/sigrw_bias_merge.log 2>&1
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh > $L/sigrw_bias_collect.log 2>&1
python3 - $BN/bias_outputs_highstat $C/output_impacts $T0 <<'PY' | tee $L/sigrw_bias_impact_summary.txt | tee -a $HB
import json, os, sys
B, I, T0 = sys.argv[1], sys.argv[2], float(sys.argv[3]); bad = 0; worst = []
for m in range(1, 31):
    p = "%s/mA_%d/merged/BiasJson/%d_gaussfit.json" % (B, m, m)
    if not (os.path.exists(p) and os.path.getmtime(p) >= T0): print("    mA%d bias MISSING" % m); bad += 1; continue
    for k, v in json.load(open(p))["fit_results"].items(): worst.append((abs(v["mean"]), m, k, v["mean"], v["mean_err"]))
    ip = "%s/%d_impacts.json" % (I, m)
    if not (os.path.exists(ip) and os.path.getmtime(ip) >= T0): print("    mA%d impact MISSING" % m); bad += 1
worst.sort(reverse=True)
print("    bias: %d truth functions; largest " % len(worst) + ", ".join("mA%d %s %+.3f+-%.3f" % w[1:] for w in worst[:4]))
over = [w for w in worst if w[0] > 0.2]
if over: print("    ABOVE 0.2: " + ", ".join("mA%d %s %+.3f" % (w[1], w[2], w[3]) for w in over)); bad += len(over)
print("VERDICT: %s" % ("PASSED" if bad == 0 else "CHECK (%d)" % bad))
PY
fi

hb "[8] DONE signal-reweight chain"
echo "$(date '+%F %H:%M') DONE chain_sigrw.sh (heartbeat sigrw_heartbeat.txt)" >> $L/PENDING_RESULTS.txt
