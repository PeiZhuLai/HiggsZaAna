#!/usr/bin/env bash
# MC-statistics bootstrap of the background-only pseudo-data closure, production vs win_mfree low-mass
# BDT (same method as pseudodata_closure/scripts/boot_mA3_closure.sh, 2026-09-27): each replica rebuilds
# the DY pseudo-data with Poisson(1) weights per MC event, refits the envelope (30 GOF toys) and reruns
# combine Significance + MultiDimFit with the given datacard (background swapped to the replica). The
# replica spread is the MC-statistics uncertainty of a single closure measurement (sd ~ 1.1 in Z at mA3),
# so production and new model are compared on the replica means, not on one value each.
#
# Usage: boot_closure.sh <prod|new> <mA> <cut> <nrep> <datacard> [npar]
# Output: /eos/home-p/pelai/HZa/sandbox_lowmass_mfree/closure_boot/<model>_mA<m>_cut<tag>/{rep_<seed>/result.txt, summary.txt}
set -u
# Fallback minimizers (not in the production combine calls): without them about 1 replica in a few has a
# failed initial fit (MultiDimFit 'failed', Significance 0). Same options for both models.
FB="--cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminFallbackAlgo Minuit2,Simplex,0:0.1 --X-rtd MINIMIZER_freezeDisassociatedParams"
MODEL=$1; MA=$2; CUT=$3; NREP=$4; DC=$5; NPAR=${6:-3}
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
CPY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
case $MODEL in
  prod) SB=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix ;;
  new)  SB=/eos/home-p/pelai/HZa/sandbox_lowmass_mfree/run3_bdt_scored ;;
  *) echo "model must be prod or new"; exit 1 ;;
esac
TAG=$(echo $CUT | tr -d '.')
OUT=/eos/home-p/pelai/HZa/sandbox_lowmass_mfree/closure_boot/${MODEL}_mA${MA}_cut${TAG}; mkdir -p $OUT
BKGDIR=$(awk '$1=="shapes" && $2=="bkg_mass" {print $4}' $DC | head -1 | xargs dirname)
[ -n "$BKGDIR" ] || { echo "cannot find the bkg_mass shapes line in $DC"; exit 1; }
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $FFIT/.. && eval `scramv1 runtime -sh`
source $FFIT/setup.sh >/dev/null 2>&1
set -u

one_rep() {
  local seed=$1 R=$OUT/rep_$1
  [ -s $R/result.txt ] && return 0
  mkdir -p $R/fit/ft
  ( cd $C/closure_scripts && env -i HOME=$HOME PATH=/usr/bin:/bin CLOSURE_FSRFIX=1 CLOSURE_SCORED_BASE=$SB \
      $CPY build_pseudodata_ws.py --mA $MA --cut $CUT --bootstrap $seed --outbase $R ) > $R/build.log 2>&1
  ( cd $FFIT/Background && ./bin/fTest_ALP_turnOn -i $R/mA_M$MA/ws/run3.root \
      --saveMultiPdf $R/fit/CMS-HGG_mva_13p6TeV_multipdf.root -D $R/fit/ft \
      --mass_ALP $MA -c 1 --isFlashgg 0 --isData 0 -f data, --datfile $R/fit/fTest.dat \
      --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 30 ) > $R/ftest.log 2>&1
  sed "s#$BKGDIR/#$R/fit/#g" "$DC" > $R/datacard.txt
  grep -q "$R/fit/CMS-HGG_mva_13p6TeV_multipdf.root" $R/datacard.txt || { echo "$seed FAIL_datacard" > $R/result.txt; return; }
  ( cd $R && text2workspace.py datacard.txt -o ws_t2w.root > t2w.log 2>&1
    combine -M Significance ws_t2w.root --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMax 2 -n B $FB > sig.log 2>&1
    combine -M MultiDimFit ws_t2w.root --algo singles --cminDefaultMinimizerStrategy 0 -m 125.38 \
        --setParameters MH=125.38 --freezeParameters MH --rMin -3 --rMax 3 -n B $FB > rhat.log 2>&1
    rm -f ws_t2w.root higgsCombine*.root )
  local Z=$(grep -E 'Significance:' $R/sig.log | sed -E 's/.*Significance: //')
  local rh=$(grep -E '^\s+r :' $R/rhat.log | awk '{print $3}')
  local N=$(grep -oE 'pseudo-data N=[0-9]+' $R/build.log | grep -oE '[0-9]+')
  grep -q "MultiDimFit failed" $R/rhat.log && rh=NA
  echo "$seed N=${N:-NA} Z=${Z:-NA} rhat=${rh:-NA}" > $R/result.txt
}
export -f one_rep; export OUT C SB CPY FFIT MA CUT DC BKGDIR FB

echo "=== $MODEL mA$MA cut $CUT: $NREP replicas, $NPAR parallel, datacard $DC ($(date '+%F %T'))"
seq 1 $NREP | xargs -P $NPAR -I{} bash -c 'one_rep {}'
cat $OUT/rep_*/result.txt | $CPY -c "
import sys, numpy as np
rows=[l.split() for l in sys.stdin if 'rhat=' in l]
r=np.array([float(x[3][5:]) for x in rows if x[3][5:] not in ('NA','')])
z=np.array([float(x[2][2:]) for x in rows if x[2][2:] not in ('NA','')])
print('$MODEL mA$MA cut $CUT: replicas ok %d/%d  rhat mean %+.4f sd %.4f (mean err %.4f)  |  Z mean %.2f sd %.2f  frac(Z>2) %.2f'
      % (r.size, len(rows), r.mean(), r.std(ddof=1), r.std(ddof=1)/np.sqrt(r.size), z.mean(), z.std(ddof=1), (z>2).mean()))" | tee $OUT/summary.txt
