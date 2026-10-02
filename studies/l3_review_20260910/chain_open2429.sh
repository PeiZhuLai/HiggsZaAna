#!/usr/bin/env bash
# mA24-30 background refit with opened step-Gaussian ranges and toy-based GOF (user OK 2026-09-27).
#
# Why: after the FSR-fix re-production the data at mA24-30 exceed 425 events, so fTest switched
# from toy GOF to the asymptotic chi2, and at mA24-29 every envelope member had its turn-on /
# width / smearing parameters railed on hard-coded per-mass windows (open-range refit dNLL up to
# 69). Changes (both already compiled into Background/bin):
#   src/PdfModelBuilder.cc  openStepGaus(): mA24-29 windows -> turnon 95-125, width 0.1-50,
#                           sigma 0.05-20 (= the bias-study fit ranges)
#   test/fTest_ALP_turnOn.cpp  GOF always by toys (only mA24-30 were on the asymptotic branch)
#
#   input : EOS root_MVAcut/data/mA_M{24..30}/ws/run3.root (unchanged)
#   output: fit_results_run3/{24..30}, root_t2w, output_combine_results, bias + impact on condor
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_fsrfix_20260925
MASSES="${MASSES:-24 25 26 27 28 29 30}"
TAG=before_open2429_$(date +%Y%m%d_%H%M)
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

echo "=============== [1] backups ($TAG) ==============="
V=$C/output_combine_results_variants; mkdir -p $V $F/Background/fit_results_variants $BN/bias_outputs_highstat_variants
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
for m in $MASSES; do cp -rp $FR/$m $F/Background/fit_results_variants/${m}_$TAG; done
echo "  done"

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
[ $F/Background/bin/fTest_ALP_turnOn -nt $F/Background/src/PdfModelBuilder.cc ] || die "fTest binary older than PdfModelBuilder.cc -- rebuild first"

echo "=============== [2] background fits (<= 6 in parallel) ==============="
T=$(date +%s)
running=0
for m in $MASSES; do
  ( cd $F/Background && bash $D/fit_bkg_single/fit_bkg_m$m.sh > $L/open2429_bkg_m$m.log 2>&1 ) &
  running=$((running+1)); [ $running -ge 6 ] && { wait -n; running=$((running-1)); }
done
wait
for m in $MASSES; do
  MP=$FR/$m/CMS-HGG_mva_13p6TeV_multipdf.root
  fresh $MP $T && [ $(stat -c %s $MP) -ge 1024 ] || die "mA$m multipdf missing/stale (open2429_bkg_m$m.log)"
  echo "  mA$m:"; sed 's/^/    /' $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt
done

echo "=============== [3] t2w + blind limits (all 30), compare ==============="
T=$(date +%s)
cd $C && sh 2_text2ws.sh > $L/open2429_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/open2429_limits.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || bad=1; done
[ $bad -eq 0 ] || die "limits stale"
python3 - $C/output_combine_results $V/output_combine_results_$TAG > $L/open2429_changed.txt <<'PY'
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
echo "  masses whose median limit moved (mA old new ratio):"; sed 's/^/    /' $L/open2429_changed.txt
CHANGED=$(awk '{print $1}' $L/open2429_changed.txt | tr '\n' ' ')
# every refit mass gets a new bias/impact even if its limit did not move (the envelope changed)
REDO=$(echo $MASSES $CHANGED | tr ' ' '\n' | sort -nu | tr '\n' ' ')

echo "=============== [4] condor: bias + impact for: $REDO ==============="
for m in $REDO; do
  cp -p $FZ/${m}_Datacard_leptons.root $FZ/${m}_Datacard_leptons.root.$TAG 2>/dev/null
  cp -p $C/root_t2w/${m}_Datacard_leptons.root $FZ/
  [ -d $BN/bias_outputs_highstat/mA_$m ] && mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG
  cp -p $C/output_impacts/${m}_impacts.json $C/output_impacts/${m}_impacts.json.$TAG 2>/dev/null
done
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $REDO 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $REDO 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "REDO_MASSES $REDO"
echo "SUBMITTED $(date '+%F %T')"
