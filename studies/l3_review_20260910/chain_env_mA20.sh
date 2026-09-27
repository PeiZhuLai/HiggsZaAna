#!/usr/bin/env bash
# mA20 envelope adjustment (user decision 2026-09-26).
# With the bias pull computed from the error on the side facing the truth (RunBiasStudy.py,
# 2026-09-26), Bern6 = +0.208 at m_a = 20 is the only function of the scan above 0.2. Per AN
# Sec 7.4 the Bernstein member falls back to the next lower order: bestBernOrder(20) 6 -> 5.
# The source edit is already in fTest_ALP_turnOn.cpp; this script rebuilds, refits m_a = 20,
# rebuilds its workspace and limit, then submits bias (10 x 1000 toys) + impact.
set -uo pipefail
STEP=${STEP:-cap5}
M=20
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
BK=$F/Background
FR=$BK/ALP_BkgModel_ReReco/fit_results_run3
C=$F/Combine
BN=$C/Checks/Bias_nominal
FZ=$C/root_t2w_fsrfix_20260925
SRC=$BK/test/fTest_ALP_turnOn.cpp
TAG=$(date +%Y%m%d_%H%M)
die() { echo "STOPPING: $*"; exit 1; }

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv

echo "=============== [1] backups (tag ${STEP}_$TAG) ==============="
[ -d $FR/$M ] && cp -rp $FR/$M $FR/${M}_bak_before_${STEP}_$TAG
cp -p $C/output_combine_results/higgsCombine$M.AsymptoticLimits.mH125.38.root $C/output_combine_results/higgsCombine$M.AsymptoticLimits.mH125.38.root.bak_before_${STEP}_$TAG
cp -p $FZ/${M}_Datacard_leptons.root $FZ/${M}_Datacard_leptons.root.bak_before_${STEP}_$TAG
# backups go OUTSIDE bias_outputs_highstat: collect_bias_results.sh globs mA_*/merged and would
# copy a variant's stale plots over the official ones (happened 2026-09-26 for mA8/mA11)
[ -d $BN/bias_outputs_highstat/mA_$M ] && mkdir -p $BN/bias_outputs_highstat_variants && mv $BN/bias_outputs_highstat/mA_$M $BN/bias_outputs_highstat_variants/mA_${M}_before_${STEP}_$TAG
echo "  done"

echo "=============== [2] source change + rebuild ==============="
grep -q 'case 20: return 5;' $SRC || die "mA20 cap not found in source"
echo "  source already capped: bestBernOrder(20) = 5"
cd $BK && make -f makefile clean > $L/env20_build.log 2>&1 && make -f makefile -B >> $L/env20_build.log 2>&1
[ -x $BK/bin/fTest_ALP_turnOn ] && [ $BK/bin/fTest_ALP_turnOn -nt $SRC ] || die "rebuild failed (env20_build.log)"
echo "  rebuilt"

echo "=============== [3] background fit m_a = $M (local) ==============="
T=$(date +%s)
bash $D/fit_bkg_single/fit_bkg_m$M.sh > $L/env20_bkg.log 2>&1
MP=$FR/$M/CMS-HGG_mva_13p6TeV_multipdf.root
[ -s $MP ] && [ $(stat -c %Y $MP) -ge $T ] && [ $(stat -c %s $MP) -ge 1024 ] || die "multipdf missing/stale/truncated"
python3 - "$MP" <<'PY' || die "multipdf check failed"
import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
f = ROOT.TFile.Open(sys.argv[1]); w = f.Get("multipdf")
mp = [p for p in w.allPdfs() if p.InheritsFrom("RooMultiPdf")][0]
names = [mp.getPdf(i).GetName() for i in range(mp.getNumPdfs())]
print("  envelope (%d): %s" % (len(names), ", ".join(names)))
sys.exit(0 if len(names) > 0 else 1)
PY
grep -A40 -i 'envelope' $FR/$M/HZAmassInde_fTest/EnvelopeResults.txt 2>/dev/null | head -25

echo "=============== [4] text2workspace + blind limit m_a = $M ==============="
cd $C
python3 RunText2Workspace.py --batch local --queue workday --mA $M > $L/env20_t2w.log 2>&1
[ $C/root_t2w/${M}_Datacard_leptons.root -nt $MP ] || die "t2w not rebuilt (env20_t2w.log)"
combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${M}_pruned_datacard_leptons.txt \
  --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind --setParameters MH=125.38 --freezeParameters MH \
  --rMax 2 -n $M > $L/env20_limit.log 2>&1
mv -f higgsCombine$M.AsymptoticLimits.mH125.38.root output_combine_results/ || die "limit failed (env20_limit.log)"
grep -E 'Expected 50' $L/env20_limit.log | sed 's/^/  /'
cp -p $C/root_t2w/${M}_Datacard_leptons.root $FZ/${M}_Datacard_leptons.root

echo "=============== [5] condor: bias (10 x 1000 toys) + impact for m_a = $M ==============="
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
