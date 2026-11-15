#!/bin/bash
# Blind Combine-review checks on the current datacards, all on the Asimov dataset (-t -1, r = 0):
#   fitdiag  MultiDimFit s+b and b-only fits with the Hesse covariance (--saveFitResult), both cards
#   rscan    MultiDimFit 1D likelihood scan of r, cross-section card, r in [0, 3 x expected median]
# Never run any of these without -t -1: on the observed data they would unblind the analysis.
# FitDiagnostics is not used: it segfaults with the --X-rtd MINIMIZER_* discrete-profiling options,
# with --robustFit 1 alone it returns no nuisance errors, and --robustHesse 1 crashes; MultiDimFit
# with --saveFitResult and the impact options gives the full covariance.
# Usage: [MASSES="1 10 30"] bash run_fitdiag_rscan_deploy.sh [fitdiag|rscan|all]   (5 in parallel). Deployed cards (11-12), output review_deploy_20261112
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/an_next_round
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
OUT=$FF/Combine/Checks/review_deploy_20261112
MODE="${1:-all}"; MASSES="${MASSES:-1 10 30}"
mkdir -p $OUT/fitdiag $OUT/rscan $W/logs
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` 2>/dev/null

fitdiag() {  # mass card: Asimov (r = 0) s+b and b-only fits with the Hesse covariance saved
  local m=$1 card=$2
  local d=$OUT/fitdiag/${card}_mA${m}
  local X="--X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
  mkdir -p $d && cd $d || return 1
  combine -M MultiDimFit $FF/Combine/root_t2w_${card}/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38 --freezeParameters MH --rMin -1 --rMax 1 --cminDefaultMinimizerStrategy 0 $X \
    --saveFitResult -n _sb > fit_sb.log 2>&1 || { echo "FAIL fit_sb $card mA$m"; return 1; }
  combine -M MultiDimFit $FF/Combine/root_t2w_${card}/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38,r=0 --freezeParameters MH,r --cminDefaultMinimizerStrategy 0 $X \
    --saveFitResult -n _b > fit_b.log 2>&1 || { echo "FAIL fit_b $card mA$m"; return 1; }
  echo "done fitdiag $card mA$m"
}
rscan() {  # mass rmax
  local m=$1 rmax=$2
  local d=$OUT/rscan/xsec_mA${m}
  mkdir -p $d && cd $d || return 1
  combine -M MultiDimFit $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38 --freezeParameters MH --algo grid --points 30 --rMin 0 --rMax $rmax \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants \
    --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 \
    -n _asimov > rscan.log 2>&1 || { echo "FAIL rscan mA$m"; return 1; }
  echo "done rscan mA$m"
}
export -f fitdiag rscan; export OUT FF

if [[ $MODE == fitdiag || $MODE == all ]]; then
  for c in xsec wilson; do for m in ${MASSES}; do echo "$m $c"; done; done | \
    xargs -P 5 -L 1 bash -c 'fitdiag $0 $1' > $W/logs/fitdiag.log 2>&1
fi
if [[ $MODE == rscan || $MODE == all ]]; then
  # r range: 3 x the expected median limit recorded at the bias submission
  awk -v M=" ${MASSES} " 'index(M, " "$1" ") {printf "%d %.4f\n", $1, 3*$2}' /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/q50_at_bias_submission.txt | \
    xargs -P 5 -L 1 bash -c 'rscan $0 $1' > $W/logs/rscan.log 2>&1
fi
date '+%F %T' > $W/fitdiag_rscan.DONE
