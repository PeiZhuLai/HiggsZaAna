#!/bin/bash
# Local blind checks for the Combine review (PLAN.md steps 1a-1d), current production workspaces.
#   asimov   b-only Asimov dataset per mass (input of the HybridNew jobs)
#   fitdiag  s+b and b-only MultiDimFit on the r=0 Asimov with the Hesse covariance (--saveFitResult), both cards
#   rscan    1D likelihood scan of r on the r=0 Asimov, xsec card, r in [0, 3 x q50]
#   statsyst fits on the r=q50 Asimov: all nuisances floating, and all constrained nuisances frozen, both cards
# Every combine call uses -t -1 (Asimov); the observed data are never used.
# FitDiagnostics is not used (segfaults with the discrete-profiling options); MultiDimFit --saveFitResult gives
# the covariance. Usage: [MASSES="1 10 20 30"] [GO=1] bash local_checks.sh [asimov|fitdiag|rscan|statsyst|all]
W=$(cd "$(dirname "$0")" && pwd)
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
OUT=$FF/Combine/Checks/review_round2
MODE="${1:-all}"; MASSES="${MASSES:-1 10 20 30}"; GO=${GO:-0}
Q50=$W/q50.txt
if [ "$GO" != 1 ]; then
  echo "DRY-RUN mode=$MODE masses=$MASSES -> $OUT/{asimov,fitdiag,rscan,statsyst}, logs $W/logs"; cat $Q50; exit 0
fi
mkdir -p $OUT/asimov $OUT/fitdiag $OUT/rscan $OUT/statsyst $W/logs
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` 2>/dev/null
X="--X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
q50() { awk -v m=$1 -v c=$2 '$1==m {print (c=="xsec") ? $2 : $3}' $Q50; }

asimov() {  # mass
  local m=$1; cd $OUT/asimov || return 1
  combine -M GenerateOnly $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --saveToys -s 123456 --setParameters MH=125.38 --freezeParameters MH -n _asimov_mA$m > gen_mA$m.log 2>&1 \
    && echo "done asimov mA$m" || echo "FAIL asimov mA$m"
}
fitdiag() {  # mass card
  local m=$1 card=$2 d=$OUT/fitdiag/${2}_mA${1}
  mkdir -p $d && cd $d || return 1
  combine -M MultiDimFit $FF/Combine/root_t2w_${card}/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38 --freezeParameters MH --rMin -1 --rMax 1 --cminDefaultMinimizerStrategy 0 $X \
    --saveFitResult -n _sb > fit_sb.log 2>&1 || { echo "FAIL fit_sb $card mA$m"; return 1; }
  combine -M MultiDimFit $FF/Combine/root_t2w_${card}/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38,r=0 --freezeParameters MH,r --cminDefaultMinimizerStrategy 0 $X \
    --saveFitResult -n _b > fit_b.log 2>&1 || { echo "FAIL fit_b $card mA$m"; return 1; }
  echo "done fitdiag $card mA$m"
}
rscan() {  # mass
  local m=$1 d=$OUT/rscan/xsec_mA${1} rmax
  rmax=$(python3 -c "print('%.5g' % (3*$(q50 $m xsec)))")
  mkdir -p $d && cd $d || return 1
  combine -M MultiDimFit $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal 0 \
    --setParameters MH=125.38 --freezeParameters MH --algo grid --points 30 --rMin 0 --rMax $rmax \
    --cminDefaultMinimizerStrategy 0 $X -n _asimov > rscan.log 2>&1 && echo "done rscan mA$m" || echo "FAIL rscan mA$m"
}
statsyst() {  # mass card
  local m=$1 card=$2 d=$OUT/statsyst/${2}_mA${1} e rr
  e=$(q50 $m $card); rr=$(python3 -c "print('%.5g' % (10*$e))")
  mkdir -p $d && cd $d || return 1
  local base="-M MultiDimFit $FF/Combine/root_t2w_${card}/${m}_Datacard_leptons.root -m 125.38 -t -1 --expectSignal $e --algo singles --robustFit 1 --cminPreScan --rMin -$rr --rMax $rr --cminDefaultMinimizerStrategy 0 $X"
  combine $base --setParameters MH=125.38 --freezeParameters MH -n _full > full.log 2>&1 || { echo "FAIL full $card mA$m"; return 1; }
  combine $base --setParameters MH=125.38 --freezeParameters MH,allConstrainedNuisances -n _stat > stat.log 2>&1 || { echo "FAIL stat $card mA$m"; return 1; }
  echo "done statsyst $card mA$m"
}
export -f asimov fitdiag rscan statsyst q50; export OUT FF X Q50

has() { [[ $MODE == all || $MODE == $1 ]]; }
has asimov   && for m in $MASSES; do echo $m; done | xargs -P 4 -L 1 bash -c 'asimov $0' >> $W/logs/local.log 2>&1
has fitdiag  && for c in xsec wilson; do for m in $MASSES; do echo "$m $c"; done; done | xargs -P 4 -L 1 bash -c 'fitdiag $0 $1' >> $W/logs/local.log 2>&1
has rscan    && for m in $MASSES; do echo $m; done | xargs -P 4 -L 1 bash -c 'rscan $0' >> $W/logs/local.log 2>&1
has statsyst && for c in xsec wilson; do for m in $MASSES; do echo "$m $c"; done; done | xargs -P 4 -L 1 bash -c 'statsyst $0 $1' >> $W/logs/local.log 2>&1
date '+%F %T' > $W/logs/local_$MODE.DONE
