#!/usr/bin/env bash
# Sandbox full chain for the win_mfree low-mass BDT at m_a = 1, 2, 3 (user request 2026-10-13).
#
# Working points (rule fixed before any limit was looked at, see closure_scan/*/summary.txt):
# closure Z < 1.0, every data-envelope GOF > 0.01, R within 2 sigma of 1, overlap share <= preselection
# and window fraction >= preselection; among passing cuts the best proxy Z.
#   mA1 0.9645   mA2 0.975   mA3 0.985        (production: 0.963 / 0.975 / 0.988)
#
# Code: a code-only copy of flashggFinalFit in ./flashggFinalFit with every path pointed at the sandbox
# (tools/commonObjects.py cwd__, root_MVAcut, scored inputs, WP JSON, S1 lnN JSON); masses restricted to
# 1-3; 3_rysn_datacard.sh (copies to ~/public/hza_datacard) disabled. The background fit uses the
# production fTest_ALP_turnOn binary with explicit sandbox input/output paths.
# Stages: [1] S1 lnN  [2] MVA cut  [3] Tree2WS  [4] background envelope  [5] signal model  [6] datacards
#         [7] blind expected limits vs production  [8] closure with the sandbox datacard  [9] bias (condor)
# Every stage is judged by its products (fresh mtimes, content checks), not exit codes. Production
# outputs are fingerprinted before and after; any change is a failure.
# Markers: logs/chain_sbx.DONE (all stages) or logs/chain_sbx.FAILED (with the reason).
set -uo pipefail
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain
S=$C/flashggFinalFit
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
SBE=/eos/home-p/pelai/HZa/sandbox_lowmass_mfree
L=$C/logs
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"
declare -A WP=([1]=0.9645 [2]=0.975 [3]=0.985)
declare -A CLOSURE_DIR=([1]=$SBE/closure_scan/new_mA1/cut_09645 [2]=$SBE/closure_scan/new_mA2/cut_0975 [3]=$SBE/closure_scan/new_mA3/cut_0985)
FR=$S/Background/ALP_BkgModel_ReReco/fit_results_run3
START=${START_STAGE:-1}
run() { [ "$START" -le "$1" ]; }
die() { echo "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/chain_sbx.FAILED; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }
rm -f $L/chain_sbx.DONE $L/chain_sbx.FAILED

fingerprint() {   # production products that must not change
  { for m in 1 2 3; do
      ls -l --time-style=+%s $F/Background/ALP_BkgModel_ReReco/fit_results_run3/$m/CMS-HGG_mva_13p6TeV_multipdf.root
      ls -l --time-style=+%s $F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt
      ls -l --time-style=+%s $F/Combine/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root
      ls -l --time-style=+%s /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M$m/run3.root
      for y in $YEARS; do ls -l --time-style=+%s /eos/home-p/pelai/HZa/root_MVAcut/sig/mA_M$m/output_$y.root; done
      ls -l --time-style=+%s $F/Signal/outdir_ele/signalFit/output/${m}_CMS-HGG_sigfit_2024_ele_Hm125.root
    done
    ls -l --time-style=+%s /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json
    ls -ld --time-style=+%s /afs/cern.ch/user/p/pelai/public/hza_datacard/mA_M1
  } 2>&1 | awk '{print $6, $7, $5}'
}
[ -s $L/prod_fingerprint.before ] || fingerprint > $L/prod_fingerprint.before
[ -e $L/chain_start.marker ] || touch $L/chain_start.marker

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }
export -f cmsenv
cmsenv
export PYTHONPATH="$S/tools:$S/Signal/tools:${PYTHONPATH:-}"
set -u
python3 -c "import commonObjects as c; assert c.cwd__ == '$S', c.cwd__; print('  [env] cwd__ =', c.cwd__)" || die "commonObjects does not point at the sandbox"

if run 1; then
echo "=============== [1/9] S1 lnN (mva_reweight) at the sandbox working points ==============="
T=$(date +%s)
( cd $C/s1 && env -i HOME=$HOME PATH=/usr/bin:/bin S1_RW_TAG=fsrfix_dyveto $ALP_PY eval_alt_reweight.py --masses 1,2,3 --out yields_sbx.csv ) > $L/s1_eval.log 2>&1 \
  || die "eval_alt_reweight failed (s1_eval.log)"
( cd $C/s1 && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY make_s1_lnn_json.py --reference nominal s1_lnn_scheme2_sbx.json yields_sbx.csv ) > $L/s1_json.log 2>&1 \
  || die "make_s1_lnn_json failed (s1_json.log)"
fresh $C/s1/s1_lnn_scheme2_sbx.json $T || die "S1 JSON not written"
cat $L/s1_json.log | sed 's/^/  /'
fi

if run 2; then
echo "=============== [2/9] MVA cut: data + signal (mA1-3, 5 eras) ==============="
T=$(date +%s)
python3 $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py > $L/mva_data.log 2>&1 || die "apply_bdt_data failed (mva_data.log)"
for y in $YEARS; do
  python3 $S/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_sig.py --samples mA_M1,mA_M2,mA_M3 --years $y > $L/mva_sig_$y.log 2>&1 &
done
wait
for m in 1 2 3; do
  fresh $SBE/root_MVAcut/data/mA_M$m/run3.root $T || die "data mA$m run3.root stale/missing"
  for y in $YEARS; do fresh $SBE/root_MVAcut/sig/mA_M$m/output_$y.root $T || die "signal mA$m $y stale/missing (mva_sig_$y.log)"; done
done
python3 - $SBE <<'PY' || die "MVA-cut trees fail a full read"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]
paths = ["%s/data/mA_M%d/run3.root" % (B + "/root_MVAcut", m) for m in (1, 2, 3)]
paths += ["%s/root_MVAcut/sig/mA_M%d/output_%s.root" % (B, m, y) for m in (1, 2, 3)
          for y in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")]
for p in paths:
    f = ROOT.TFile.Open(p); n = 0; ntree = 0
    def walk(d):
        global n, ntree
        for k in d.GetListOfKeys():
            o = k.ReadObj()
            if o.InheritsFrom("TDirectory"): walk(o)
            elif o.InheritsFrom("TTree"):
                ntree += 1
                for i in range(o.GetEntries()):
                    if o.GetEntry(i) <= 0: print("FAIL", p, o.GetName(), i); sys.exit(1)
                n += o.GetEntries()
    walk(f)
    print("  %-90s trees=%3d entries=%d" % (p.replace(B, "$SBE"), ntree, n))
PY
fi

if run 3; then
echo "=============== [3/9] Tree2WS ==============="
T=$(date +%s)
( cd $S/Trees2WS && bash run_tree2ws.sh ) > $L/tree2ws.log 2>&1
for m in 1 2 3; do
  fresh $SBE/root_MVAcut/data/mA_M$m/ws/run3.root $T || die "data ws mA$m stale/missing (tree2ws.log)"
  n=$(find $SBE/root_MVAcut/sig/mA_M$m/ws_Tree2WS -name '*.root' -newermt "@$T" 2>/dev/null | wc -l)
  [ "$n" -ge 10 ] || die "signal ws mA$m: only $n fresh files"
done
# the data workspaces must hold exactly the events the closure scan fitted at the same cut
python3 - $SBE ${CLOSURE_DIR[1]} ${CLOSURE_DIR[2]} ${CLOSURE_DIR[3]} <<'PY' || die "data workspace counts differ from the closure-scan data"
import sys, re, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
B = sys.argv[1]; bad = 0
for m, cd in zip((1, 2, 3), sys.argv[2:]):
    f = ROOT.TFile.Open("%s/root_MVAcut/data/mA_M%d/ws/run3.root" % (B, m))
    ws = f.Get("CMS_hza_workspace"); d = ws.data("Data_13p6TeV")
    n = d.sumEntries()
    want = int(re.search(r"data N=(\d+)", open(cd + "/data.build.log").read()).group(1))
    print("  mA%d data ws entries %d, closure-scan data N %d %s" % (m, n, want, "OK" if n == want else "MISMATCH"))
    bad |= (n != want)
sys.exit(bad)
PY
fi

if run 4; then
echo "=============== [4/9] background envelopes (production binary, sandbox paths, 3 in parallel) ==============="
T=$(date +%s)
bkg_one() {
  local m=$1 P=$FR/$1
  mkdir -p $P/HZAmassInde_fTest $FR/AllFitResults_blind $FR/AllFitResults_unblind
  ( cd $F/Background && ./bin/fTest_ALP_turnOn -i $SBE/root_MVAcut/data/mA_M$m/ws/run3.root \
      --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root -D $P/HZAmassInde_fTest \
      --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $P/fTest.dat \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  for b in blind unblind; do
    extra=""; [ $b = unblind ] && extra="--unblind"
    ( cd $F/Background && ./bin/makeBkgPlots_ALP -b $P/CMS-HGG_mva_13p6TeV_multipdf.root -d $P/BkgPlots_$b \
        --total_OutDir $FR/AllFitResults_$b -o $P/BkgPlots_$b.root --sqrts 13p6TeV --isMultiPdf --useBinnedData \
        --massStep 2.5 --mhVal 125.0 --maVal $m --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
        --intLumi 172.13 -c 0 --isFlashgg 0 --doBands $extra ) >> $P/makeBkgPlots.stdout.log 2>&1
  done
}
for m in 1 2 3; do bkg_one $m & done; wait
for m in 1 2 3; do
  MP=$FR/$m/CMS-HGG_mva_13p6TeV_multipdf.root
  fresh $MP $T && [ $(stat -c %s $MP) -ge 1024 ] || die "mA$m multipdf missing/stale/truncated"
  echo "  mA$m: $(grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')"
  python3 - $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt <<'PY' || die "mA$m: a data-envelope member has GOF <= 0.01"
import re, sys
g = [float(x) for x in re.findall(r"gof : ([0-9.e-]+)", open(sys.argv[1]).read())]
sys.exit(0 if g and min(g) > 0.01 else 1)
PY
done
fi

if run 5; then
echo "=============== [5/9] signal model (fTest -> photon syst -> signalFit -> plotter), one mass per process ==============="
T=$(date +%s)
mkdir -p $C/sig_single
sig_one() {
  local m=$1
  for s in 1_runjob_sig_fTest 2_runjob_sig_calcPhotonSyst 3_runjob_sig_signalFit 4_runjob_sig_RunPlotter; do
    sed -e "s/^mAs=( 1 2 3 )$/mAs=( $m )/" $S/shellScripts/sig_sys/$s.sh > $C/sig_single/${s}_m$m.sh
    grep -q "^mAs=( $m )$" $C/sig_single/${s}_m$m.sh || { echo "copy of $s failed"; return 1; }
    echo "  [$(date '+%H:%M:%S')] mA$m $s"
    bash $C/sig_single/${s}_m$m.sh > $L/sig_${s}_m$m.log 2>&1
  done
}
for m in 1 2 3; do sig_one $m & done; wait
for m in 1 2 3; do for y in $YEARS; do for ch in ele mu; do
  fresh $S/Signal/outdir_$ch/signalFit/output/${m}_CMS-HGG_sigfit_${y}_${ch}_Hm125.root $T || die "signal fit mA$m $y $ch stale/missing (sig_*_m$m.log)"
done; done; done
echo "  signal fits: 30/30 fresh"
fi

if run 6; then
echo "=============== [6/9] datacards (makeYields + makeDatacard; rysn disabled) ==============="
T=$(date +%s)
( cd $S/Datacard && sh 1_runjob_gen_datacard_makeYields.sh ) > $L/dc_1.log 2>&1
( cd $S/Datacard && sh 2_runjob_gen_datacard_makeDatacard.sh ) > $L/dc_2.log 2>&1
for m in 1 2 3; do
  DC=$S/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt
  fresh $DC $T || die "datacard mA$m stale/missing (dc_*.log)"
  # every shape must come from the sandbox
  if awk '$1=="shapes"{for(i=1;i<=NF;i++) if($i ~ /\.root$/) print $i}' $DC | grep -v "^$S/" | grep -q .; then
    awk '$1=="shapes"{for(i=1;i<=NF;i++) if($i ~ /\.root$/) print $i}' $DC | grep -v "^$S/" | sort -u | head
    die "datacard mA$m references shapes outside the sandbox"
  fi
  echo "  mA$m datacard OK ($(grep -c '^shapes' $DC) shapes lines, all in the sandbox); mva_reweight: $(grep -E '^CMS_hza_mva_reweight' $DC | awk '{print $2, $3}')"
done
fi

if run 7; then
echo "=============== [7/9] text2workspace + blind expected limits, compared with production ==============="
T=$(date +%s)
mkdir -p $S/Combine/root_t2w $S/Combine/output_combine_results
for m in 1 2 3; do
  ( cd $S/Combine && text2workspace.py $S/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt \
      -o root_t2w/${m}_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > $L/t2w_m$m.log 2>&1
    combine -M AsymptoticLimits root_t2w/${m}_Datacard_leptons.root --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind \
      --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n $m > $L/limit_m$m.log 2>&1
    mv higgsCombine$m.AsymptoticLimits.mH125.38.root output_combine_results/ ) &
done
wait
for m in 1 2 3; do fresh $S/Combine/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || die "limit mA$m missing (limit_m$m.log)"; done
python3 - $S/Combine/output_combine_results $F/Combine/output_combine_results <<'PY' | tee $C/limits_vs_production.txt
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def q(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); out = {}
    for e in f.Get("limit"): out[round(e.quantileExpected, 3)] = e.limit
    f.Close(); return out
print("expected limit on r (blind Asimov): production vs sandbox win_mfree")
print("  mA   prod -1s   prod med   prod +1s |  sbx -1s    sbx med    sbx +1s | med sbx/prod")
for m in (1, 2, 3):
    p, s = q(sys.argv[2], m), q(sys.argv[1], m)
    print("  %d  %9.4f  %9.4f  %9.4f | %8.4f  %9.4f  %9.4f | %8.3f" % (m, p[0.16], p[0.5], p[0.84], s[0.16], s[0.5], s[0.84], s[0.5] / p[0.5]))
PY
fi

if run 8; then
echo "=============== [8/9] closure with the sandbox datacard (new signal model) ==============="
for m in 1 2 3; do
  CD=${CLOSURE_DIR[$m]}; O=$C/closure_final/mA$m; mkdir -p $O
  sed "s#$FR/$m/#$CD/ps/fit/#g" $S/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt > $O/datacard.txt
  grep -q "$CD/ps/fit/CMS-HGG_mva_13p6TeV_multipdf.root" $O/datacard.txt || die "closure datacard mA$m: background path not swapped"
  ( cd $O && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n PS$m > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n PS$m > rhat.log 2>&1 ) &
done
wait
for m in 1 2 3; do
  O=$C/closure_final/mA$m
  echo "  mA$m closure Z=$(grep -E 'Significance:' $O/sig.log | sed -E 's/.*Significance: //')  r_hat=$(grep -E '^\s+r :' $O/rhat.log | sed -E 's/^\s+r :\s+//; s/ \(68%\)//')"
done | tee $C/closure_final/summary.txt
fi

if run 9; then
echo "=============== [9/9] bias (high-stat, condor, 10 x 1000 toys per mass, batch Bias_HighStat_sbx) ==============="
ROOT_DATACARD_PATH=$S/Combine/root_t2w NCHUNKS=10 bash $S/shellScripts/bias/Condor/subjob_bias_highstat.sh 1 2 3 2>&1 | tail -5 | sed 's/^/  /'
fi

fingerprint > $L/prod_fingerprint.after
diff $L/prod_fingerprint.before $L/prod_fingerprint.after > $L/prod_fingerprint.diff && echo "  production products unchanged" \
  || die "production products changed: see prod_fingerprint.diff"
# anything written into the production areas since the chain started is a failure
timeout 900 find $F /eos/home-p/pelai/HZa/root_MVAcut -type f -newer $L/chain_start.marker > $L/prod_new_files.txt 2>/dev/null
[ -s $L/prod_new_files.txt ] && { head $L/prod_new_files.txt; die "files written into production areas (prod_new_files.txt)"; }
echo "  no file in production areas newer than the chain start"
date '+%F %T' > $L/chain_sbx.DONE
echo "=============== SANDBOX CHAIN DONE $(date '+%F %T') ==============="
