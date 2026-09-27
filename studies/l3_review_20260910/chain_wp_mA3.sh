#!/usr/bin/env bash
# mA3 working point 0.988 -> 0.940 (user decision 2026-09-26).
#
# The low-mass working points follow the AN's R ~ 1 criterion, with the tolerance used for mA2
# on 2026-07-26 (R 1.25 -> ~1.1). After the FSR-fix retraining the pinned mA3 cut 0.988 sits at
# R = 1.21. R rises monotonically above cut ~0.90, so the tightest cut with R <= 1.10 is 0.940
# (R 1.100, Asimov Z 26.9 vs 28.8 at 0.988). The strict R = 1 crossing is at 0.721 (Z 22.1).
#
# Only mA3 changes at the MVA-cut / workspace / background / signal level (interpolated points
# start at mA11 and take their cut from the nearest anchor, never mA3). The efficiency JSON is
# regenerated and ALL datacards + limits are redone, because the quadratic efficiency curve
# through the anchors could move the interpolated points; the chain reports which limits moved.
#
#   input : Plot/output/MVAcut_points_run3.json, run3_bdt_scored_fsrfix, EOS root_MVAcut
#   output: EOS root_MVAcut/{data,sig}/mA_M3, fit_results_run3/3, Signal outdir 3_*, datacards,
#           Combine/output_combine_results, root_t2w(+ frozen copy), bias + impact on condor
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
PLOT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_fsrfix_20260925
EOSMVA=/eos/home-p/pelai/HZa/root_MVAcut
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
NEWCUT=0.940
TAG=before_wp0940_$(date +%Y%m%d_%H%M)
S1=$D/wp3_single; mkdir -p $S1
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

echo "=============== [1] backups ($TAG) ==============="
V=$C/output_combine_results_variants; mkdir -p $V
# 1_makeLimits.sh does rm -fr output_combine_results: move every non-standard file out first
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
cp -p $PLOT/output/MVAcut_points_run3.json $PLOT/output/MVAcut_points_run3.json.$TAG
cp -p $PLOT/scripts/scan_score_R_significance.py $PLOT/scripts/scan_score_R_significance.py.$TAG
for d in data sig; do [ -d $EOSMVA/$d/mA_M3 ] && cp -rp $EOSMVA/$d/mA_M3 $EOSMVA/${d}_variants_mA_M3_$TAG; done
# variants live OUTSIDE ALP_BkgModel_ReReco: sync_figures.sh copies that whole tree into the AN
mkdir -p $F/Background/fit_results_variants && cp -rp $FR/3 $F/Background/fit_results_variants/3_$TAG
mkdir -p $F/Signal/sigfit_variants_$TAG && cp -p $F/Signal/outdir_{ele,mu}/signalFit/output/3_* $F/Signal/sigfit_variants_$TAG/ 2>/dev/null
cp -rp $F/Datacard/output_Datacard_leptons $F/Datacard/output_Datacard_leptons_$TAG
cp -p $FZ/3_Datacard_leptons.root $FZ/3_Datacard_leptons.root.$TAG
mkdir -p $BN/bias_outputs_highstat_variants && [ -d $BN/bias_outputs_highstat/mA_3 ] && mv $BN/bias_outputs_highstat/mA_3 $BN/bias_outputs_highstat_variants/mA_3_$TAG
cp -p $C/output_impacts/3_impacts.json $C/output_impacts/3_impacts.json.$TAG 2>/dev/null
echo "  done"

echo "=============== [2] working point: FIXED_WP[3] -> $NEWCUT, rewrite JSON ==============="
python3 - $PLOT/scripts/scan_score_R_significance.py $NEWCUT <<'PY' || exit 1
import sys
p, cut = sys.argv[1], sys.argv[2]; s = open(p).read()
old = "FIXED_WP = {2: 0.975, 3: 0.988}"
new = ("FIXED_WP = {2: 0.975, 3: %s}   # mA3 0.988 -> %s (2026-09-26): after the FSR-fix retraining\n"
       "                                  # 0.988 sat at R = 1.21; %s is the tightest cut with R <= 1.10") % (cut, cut, cut)
assert s.count(old) == 1, "FIXED_WP line not found"
open(p, "w").write(s.replace(old, new)); print("  FIXED_WP updated")
PY
cd $PLOT && PYTHONPATH=$PLOT/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts $ALP_PY scripts/scan_score_R_significance.py --write-json > $L/wp3_scan.log 2>&1
grep -E '^ +[0-9]+ +0\.' $L/wp3_scan.log
python3 - $PLOT/output/MVAcut_points_run3.json $NEWCUT <<'PY' || die "JSON not updated"
import json, sys
r = {e["mA"]: e["MVAcut"] for e in json.load(open(sys.argv[1]))["results"]}
print("  JSON mA1-4:", {m: r[m] for m in (1, 2, 3, 4)})
sys.exit(0 if abs(r[3] - float(sys.argv[2])) < 1e-6 and r[2] == 0.975 else 1)
PY

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"

echo "=============== [3] MVA cut mA3: data (full-read verified) + signal 5 eras ==============="
bash $D/regen_mvacut_data.sh 3 > $L/wp3_mva_data.log 2>&1 || die "data MVA cut failed (wp3_mva_data.log)"
tail -2 $L/wp3_mva_data.log
T=$(date +%s)
for y in $YEARS; do
  python3 $F/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M3 --years $y > $L/wp3_mva_sig_$y.log 2>&1 &
done
wait
for y in $YEARS; do fresh $EOSMVA/sig/mA_M3/output_$y.root $T || die "signal MVA cut mA3 $y stale/missing"; done
python3 - $EOSMVA $YEARS <<'PY' || die "signal MVA-cut trees fail a full read"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B, years = sys.argv[1], sys.argv[2:]
for y in years:
    f = ROOT.TFile.Open("%s/sig/mA_M3/output_%s.root" % (B, y)); n = 0
    for k in f.GetListOfKeys():
        o = k.ReadObj()
        if o.InheritsFrom("TDirectory"):
            for kk in o.GetListOfKeys():
                t = kk.ReadObj()
                if t.InheritsFrom("TTree"):
                    for i in range(t.GetEntries()):
                        if t.GetEntry(i) <= 0: print("FAIL", y, t.GetName(), i); sys.exit(1)
                    n += t.GetEntries()
    print("  sig mA3 %-13s full read OK (%d entries over all trees)" % (y, n))
PY

echo "=============== [4] Tree2WS mA3 ==============="
sed -e 's/^mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)$/mAs_sig=(3)/' \
    -e 's/^mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)$/mAs_data=(3)/' \
    $F/Trees2WS/run_tree2ws.sh > $S1/run_tree2ws_m3.sh
grep -q '^mAs_sig=(3)$' $S1/run_tree2ws_m3.sh && grep -q '^mAs_data=(3)$' $S1/run_tree2ws_m3.sh || die "tree2ws single-mass copy failed"
T=$(date +%s)
cd $F/Trees2WS && bash $S1/run_tree2ws_m3.sh > $L/wp3_tree2ws.log 2>&1
fresh $EOSMVA/data/mA_M3/ws/run3.root $T || die "data ws mA3 stale"
n=$(find $EOSMVA/sig/mA_M3/ws_Tree2WS -name '*.root' -newermt "@$T" | wc -l); [ $n -ge 10 ] || die "only $n fresh signal ws for mA3"
echo "  data ws + $n signal ws fresh"

echo "=============== [5] background fit mA3 (local) ==============="
T=$(date +%s)
bash $D/fit_bkg_single/fit_bkg_m3.sh > $L/wp3_bkg.log 2>&1
MP=$FR/3/CMS-HGG_mva_13p6TeV_multipdf.root
fresh $MP $T && [ $(stat -c %s $MP) -ge 1024 ] || die "mA3 multipdf missing/stale/truncated"
cat $FR/3/HZAmassInde_fTest/EnvelopeResults.txt | sed 's/^/  /'

echo "=============== [6] signal model mA3 (fTest -> syst -> fit -> plotter -> effSigma) ==============="
T=$(date +%s)
for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
  sed -e 's/^mAs=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )$/mAs=( 3 )/' -e 's/^mAs_lowMA=( .* )$/mAs_lowMA=( )/' \
      $F/shellScripts/sig_sys/$s.sh > $S1/$s.sh
  grep -q '^mAs=( 3 )$' $S1/$s.sh || die "single-mass copy of $s failed"
  echo "  [$(date '+%H:%M:%S')] $s"
  bash $S1/$s.sh > $L/wp3_sig_$s.log 2>&1
done
bash $F/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh > $L/wp3_sig_5.log 2>&1
for y in $YEARS; do for ch in ele mu; do
  fresh $F/Signal/outdir_$ch/signalFit/output/3_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || die "signal fit mA3 $y $ch stale"
done; done
echo "  signal fits mA3: 10/10 fresh"

echo "=============== [7] efficiency JSON + all datacards ==============="
( cd $PLOT && env -u PYTHONPATH -u PYTHONHOME PYTHONPATH=$PLOT/lib $ALP_PY scripts/signal_eff_sumw.py > $L/wp3_signal_eff_sumw.log 2>&1 ) || die "signal_eff_sumw failed"
T=$(date +%s)
cd $F/Datacard
sh 1_runjob_gen_datacard_makeYields.sh  > $L/wp3_dc_1.log 2>&1
sh 2_runjob_gen_datacard_makeDatacard.sh > $L/wp3_dc_2.log 2>&1
sh 3_rysn_datacard.sh                    > $L/wp3_dc_3.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt $T || bad=1; done
[ $bad -eq 0 ] || die "datacards stale"
for m in $(seq 1 30); do
  n=$(grep -c rateParam $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt)
  case " 1 2 3 4 5 6 7 8 9 10 15 20 25 30 " in *" $m "*) want=0;; *) want=10;; esac
  [ "$n" -eq "$want" ] || die "mA$m has $n rateParam lines (want $want)"
done
echo "  datacards 30/30 fresh, rateParams OK"

echo "=============== [8] t2w + blind limits (all 30), compare ==============="
T=$(date +%s)
cd $C && sh 2_text2ws.sh > $L/wp3_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/wp3_limits.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || bad=1; done
[ $bad -eq 0 ] || die "limits stale"
python3 - $C/output_combine_results $V/output_combine_results_$TAG > $L/wp3_changed.txt <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for m in range(1, 31):
    n, o = med(sys.argv[1], m), med(sys.argv[2], m)
    if abs(n / o - 1) > 1e-4: print(m, "%.5f %.5f %.4f" % (o, n, n / o))
PY
echo "  masses whose median limit moved (mA old new ratio):"; sed 's/^/    /' $L/wp3_changed.txt
CHANGED=$(awk '{print $1}' $L/wp3_changed.txt | tr '\n' ' ')
echo "$CHANGED" | grep -qw 3 || die "mA3 limit did not change -- the new working point did not propagate"

echo "=============== [9] condor: bias + impact for changed masses: $CHANGED ==============="
for m in $CHANGED; do
  cp -p $FZ/${m}_Datacard_leptons.root $FZ/${m}_Datacard_leptons.root.$TAG 2>/dev/null
  cp -p $C/root_t2w/${m}_Datacard_leptons.root $FZ/
  [ -d $BN/bias_outputs_highstat/mA_$m ] && mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG
done
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $CHANGED 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $CHANGED 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
