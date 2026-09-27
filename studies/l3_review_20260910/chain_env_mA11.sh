#!/usr/bin/env bash
# mA11 envelope adjustment, following AN Sec 7.4 (user decision 2026-09-25).
#
# After the FSR fix the 10000-toy bias study gave Bern3 = -0.288 at m_a = 11 -- the only
# function of that point above the 0.2 threshold. The AN procedure: if the Bernstein member is
# the one above threshold, fall back to a lower order; if no lower order is available (or it
# still fails), drop the family, as done at m_a = 2, 4, 27. At m_a = 11 the low-order envelope
# allowed Bern1..Bern3, so the next step is Bern <= 2.
#
#   STEP=cap2  (default) cap Bernstein at order 2 -> refit, t2w, limit, bias + impact on condor
#   STEP=drop            remove Bernstein at m_a = 11 (use if Bern2 also fails)
#
# Only m_a = 11 is refit. The rebuilt binary is shared, but the source change only touches the
# m_a = 11 branch. Bias reads the frozen workspace copy (root_t2w_fsrfix_20260925), because the
# impact job re-runs text2workspace into Combine/root_t2w/11.
set -uo pipefail
STEP=${STEP:-cap2}
M=11
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
if [ "$STEP" = cap2 ]; then
  python3 - "$SRC" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
old = '    if (funcType=="Bernstein") return 4;   // Bern1..Bern3\n'
new = ('    // [PZ 2026-09-25] after the FSR fix Bern3 = -0.288 at m_a = 11 (10000 toys), the only\n'
       '    // function above threshold; per AN Sec 7.4 fall back to the next lower order.\n'
       '    if (funcType=="Bernstein") return (mass_ALP==11) ? 3 : 4;   // m_a=11: Bern1..Bern2\n')
assert s.count(old) == 1, "envMaxOrder Bernstein line not found"
open(p, "w").write(s.replace(old, new)); print("  envMaxOrder: m_a=11 Bernstein capped at order 2")
PY
elif [ "$STEP" = drop ]; then
  python3 - "$SRC" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
old = '    case 10: return 5; case 11: return 5; case 12: return 5; case 13: return 4;\n'
new = ('    // [PZ 2026-09-25] m_a = 11 drops Bernstein: after the FSR fix Bern3 and Bern2 both\n'
       '    // exceed the 0.2 bias threshold at 10000 toys (AN Sec 7.4).\n'
       '    case 10: return 5; case 11: return 0; case 12: return 5; case 13: return 4;\n')
assert s.count(old) == 1, "bestBernOrder case 11 not found"
open(p, "w").write(s.replace(old, new)); print("  bestBernOrder: m_a=11 -> 0 (Bernstein dropped)")
PY
else
  die "unknown STEP=$STEP"
fi
[ $? -eq 0 ] || die "source edit failed"
cd $BK && make -f makefile clean > $L/env11_build.log 2>&1 && make -f makefile -B >> $L/env11_build.log 2>&1
[ -x $BK/bin/fTest_ALP_turnOn ] && [ $BK/bin/fTest_ALP_turnOn -nt $SRC ] || die "rebuild failed (env11_build.log)"
echo "  rebuilt"

echo "=============== [3] background fit m_a = $M (local) ==============="
T=$(date +%s)
bash $D/fit_bkg_single/fit_bkg_m$M.sh > $L/env11_bkg.log 2>&1
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
python3 RunText2Workspace.py --batch local --queue workday --mA $M > $L/env11_t2w.log 2>&1
[ $C/root_t2w/${M}_Datacard_leptons.root -nt $MP ] || die "t2w not rebuilt (env11_t2w.log)"
combine -M AsymptoticLimits $F/Datacard/output_Datacard_leptons/${M}_pruned_datacard_leptons.txt \
  --cminDefaultMinimizerStrategy 0 -m 125.38 --run blind --setParameters MH=125.38 --freezeParameters MH \
  --rMax 2 -n $M > $L/env11_limit.log 2>&1
mv -f higgsCombine$M.AsymptoticLimits.mH125.38.root output_combine_results/ || die "limit failed (env11_limit.log)"
grep -E 'Expected 50' $L/env11_limit.log | sed 's/^/  /'
cp -p $C/root_t2w/${M}_Datacard_leptons.root $FZ/${M}_Datacard_leptons.root

echo "=============== [5] condor: bias (10 x 1000 toys) + impact for m_a = $M ==============="
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $M 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
