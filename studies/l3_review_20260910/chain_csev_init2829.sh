#!/usr/bin/env bash
# (a) mA28/29 background refit with the Exp-family starting values moved to where the fits end
#     (turn-on 107, width 4, sigma 9, slope -0.07; slopes down to -0.20) -- tested 2026-09-27:
#     mA28 Exp1 now converges (GOF 0.88, no family fallback), mA29 Exp1 converges (GOF 0.98).
# (b) CSEV nuisance CMS_hza_photon_csev (per era) added to Datacard/systematics_HToZa.py: the
#     weight was applied centrally but its uncertainty was missing from every datacard.
# Then all 30 datacards, t2w, blind limits; bias for mA28/29, expected impacts for all 30.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
CMSSW_SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
F=$CMSSW_SRC/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
FZ=$C/root_t2w_fsrfix_20260925
TAG=before_csev2829_$(date +%Y%m%d_%H%M)
die() { echo "STOPPING: $*"; exit 1; }
fresh() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -ge "$2" ]; }

echo "=============== [1] backups ($TAG) ==============="
V=$C/output_combine_results_variants; mkdir -p $V $F/Background/fit_results_variants
find $C/output_combine_results -maxdepth 1 -type f ! -regex '.*/higgsCombine[0-9]+\.AsymptoticLimits\.mH125\.38\.root' -exec mv {} $V/ \;
cp -rp $C/output_combine_results $V/output_combine_results_$TAG
for m in 28 29; do cp -rp $FR/$m $F/Background/fit_results_variants/${m}_$TAG; done
cp -rp $F/Datacard/output_Datacard_leptons $F/Datacard/output_Datacard_leptons_$TAG
echo "  done"

set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; set -u
cmsenv() { eval "$(cd $CMSSW_SRC && scramv1 runtime -sh)"; }; export -f cmsenv
export PYTHONPATH="${PYTHONPATH:-}:$F/tools:$F/Signal/tools"
[ $F/Background/bin/fTest_ALP_turnOn -nt $F/Background/src/PdfModelBuilder.cc ] || die "fTest binary older than PdfModelBuilder.cc"

echo "=============== [2] background fits mA28, mA29 ==============="
T=$(date +%s)
for m in 28 29; do ( cd $F/Background && bash $D/fit_bkg_single/fit_bkg_m$m.sh > $L/csev2829_bkg_m$m.log 2>&1 ) & done
wait
for m in 28 29; do
  fresh $FR/$m/CMS-HGG_mva_13p6TeV_multipdf.root $T || die "mA$m multipdf stale"
  echo "  mA$m:"; sed 's/^/    /' $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt
done

echo "=============== [3] all datacards (with CMS_hza_photon_csev) ==============="
T=$(date +%s)
cd $F/Datacard
sh 1_runjob_gen_datacard_makeYields.sh  > $L/csev2829_dc_1.log 2>&1
sh 2_runjob_gen_datacard_makeDatacard.sh > $L/csev2829_dc_2.log 2>&1
sh 3_rysn_datacard.sh                    > $L/csev2829_dc_3.log 2>&1
for m in $(seq 1 30); do
  dc=$F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons.txt
  fresh $dc $T || die "datacard mA$m stale"
  n=$(grep -c "^CMS_hza_photon_csev" $dc); [ "$n" -ge 1 ] || die "mA$m has no CMS_hza_photon_csev line"
  r=$(grep -c rateParam $dc)
  case " 1 2 3 4 5 6 7 8 9 10 15 20 25 30 " in *" $m "*) want=0;; *) want=10;; esac
  [ "$r" -eq "$want" ] || die "mA$m has $r rateParam lines (want $want)"
done
grep "^CMS_hza_photon_csev" $F/Datacard/output_Datacard_leptons/5_pruned_datacard_leptons.txt | cut -c1-160 | sed 's/^/    /'
echo "  datacards 30/30 fresh, CSEV present, rateParams OK"

echo "=============== [4] t2w + blind limits (all 30), compare ==============="
T=$(date +%s)
cd $C && sh 2_text2ws.sh > $L/csev2829_t2w.log 2>&1 && sh 1_makeLimits.sh > $L/csev2829_limits.log 2>&1
bad=0; for m in $(seq 1 30); do fresh $C/output_combine_results/higgsCombine$m.AsymptoticLimits.mH125.38.root $T || bad=1; done
[ $bad -eq 0 ] || die "limits stale"
python3 - $C/output_combine_results $V/output_combine_results_$TAG > $L/csev2829_changed.txt <<'PY'
import sys, ROOT
ROOT.gROOT.SetBatch(True)
def med(d, m):
    f = ROOT.TFile.Open("%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (d, m)); t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for m in range(1, 31):
    n, o = med(sys.argv[1], m), med(sys.argv[2], m)
    print(m, "%.5f %.5f %.4f" % (o, n, n / o))
PY
echo "  mA old new ratio:"; sed 's/^/    /' $L/csev2829_changed.txt

echo "=============== [5] condor: bias mA28 mA29, expected impacts all 30 ==============="
mkdir -p $BN/bias_outputs_highstat_variants
# only 28/29 here: the bias jobs of 24-27,30 still read $FZ; the others are refreshed after they finish
for m in 28 29; do
  cp -p $FZ/${m}_Datacard_leptons.root $FZ/${m}_Datacard_leptons.root.$TAG 2>/dev/null
  cp -p $C/root_t2w/${m}_Datacard_leptons.root $FZ/
done
for m in 28 29; do [ -d $BN/bias_outputs_highstat/mA_$m ] && mv $BN/bias_outputs_highstat/mA_$m $BN/bias_outputs_highstat_variants/mA_${m}_$TAG; done
ROOT_DATACARD_PATH=$FZ NCHUNKS=10 bash $F/shellScripts/bias/Condor/subjob_bias_highstat.sh 28 29 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  bias: /'
bash $F/shellScripts/impact/Condor/subjob_expectedImpact.sh $(seq 1 30) 2>&1 | grep -E 'submitted|ERROR' | sed 's/^/  impact: /'
echo "SUBMITTED $(date '+%F %T')"
