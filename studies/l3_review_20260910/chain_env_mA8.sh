#!/usr/bin/env bash
# mA8 envelope adjustment (user decision 2026-09-25: remove Pow1 first).
#
# After the FSR fix the 10000-toy bias study gave Pow1 = -0.223 and Lau1 = -0.204 at m_a = 8.
# The envelope is already Exp1/Pow1/Lau1 (lowest orders, no Bernstein), so the AN Sec 7.4
# fallback-to-lower-order rule has nothing to act on. Following the AN precedent at m_a = 2
# (removing the member that absorbs signal also reduced the others' pulls), the worst member,
# Pow1, is removed from the envelope and the bias study is repeated. If Lau1 is still above
# 0.2 afterwards, report back -- do not remove it too (that would leave a one-function envelope).
#
# Only m_a = 8 is refit. Bias reads the frozen workspace copy (root_t2w_fsrfix_20260925).
set -uo pipefail
STEP=${STEP:-nopow}
M=8
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
cp -p $SRC $SRC.bak_before_${STEP}_$TAG
echo "  done"

echo "=============== [2] source change + rebuild ==============="
if [ "$STEP" = nopow ]; then
  python3 - "$SRC" <<'PY2'
import sys
p = sys.argv[1]; s = open(p).read()
old = '\tfunctionClasses.push_back("PowerLaw");\n'
new = ('\t// [PZ 2026-09-25] m_a = 8: after the FSR fix Pow1 = -0.223 and Lau1 = -0.204 (10000 toys),\n'
       '\t// with every member already at the lowest order. Pow1, the worst, is removed (AN Sec 7.4).\n'
       '\tif (mass_ALP != 8) functionClasses.push_back("PowerLaw");\n')
assert s.count(old) == 1, "PowerLaw push_back not found"
open(p, "w").write(s.replace(old, new)); print("  functionClasses: PowerLaw removed at m_a=8")
PY2
else
  die "unknown STEP=$STEP"
fi
[ $? -eq 0 ] || die "source edit failed"
cd $BK && make -f makefile clean > $L/env8_build.log 2>&1 && make -f makefile -B >> $L/env8_build.log 2>&1
[ -x $BK/bin/fTest_ALP_turnOn ] && [ $BK/bin/fTest_ALP_turnOn -nt $SRC ] || die "rebuild failed (env8_build.log)"
echo "  rebuilt"

echo "=============== [3] background fit m_a = $M (local) ==============="
T=$(date +%s)
bash $D/fit_bkg_single/fit_bkg_m$M.sh > $L/env8_bkg.log 2>&1
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
python3 RunText2Workspace.py --batch local --queue workday --mA $M > $L/env8_t2w.log 2>&1
[ $C/root_t2w/${M}_Datacard_leptons.root -nt $MP ] || die "t2w not rebuilt (env8_t2w.log)"
combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${M}_pruned_datacard_leptons.txt \
  --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind --setParameters MH=125.38 --freezeParameters MH \
  --rMax 2 -n $M > $L/env8_limit.log 2>&1
mv -f higgsCombine$M.AsymptoticLimits.mH125.38.root output_combine_results/ || die "limit failed (env8_limit.log)"
grep -E 'Expected 50' $L/env8_limit.log | sed 's/^/  /'
cp -p $C/root_t2w/${M}_Datacard_leptons.root $FZ/${M}_Datacard_leptons.root

echo "=============== [5] condor: bias (10 x 1000 toys) + impact for m_a = $M ==============="
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
