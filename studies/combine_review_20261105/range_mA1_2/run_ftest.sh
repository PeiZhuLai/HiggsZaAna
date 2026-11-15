#!/usr/bin/env bash
# fTest of mA1 and mA2 with the staging build (turn-on 95-113 GeV, width and smearing opened, several starting points),
# production command line of shellScripts/bkg/fit_bkg.sh (mhLow 100 for mA1, 95 otherwise; 100 GOF toys).
# Input : /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M<m>/ws/run3.root
# Output: fit_results/<m>/{CMS-HGG_mva_13p6TeV_multipdf.root,HZAmassInde_fTest/,ftest.stdout.log,.ftest_rc}; logs/ftest.DONE
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd); S=$D/Background
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
MV=/eos/home-p/pelai/HZa/root_MVAcut/data
rm -f $D/logs/ftest.DONE
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
one() {
  local m=$1 lo=95 P=$D/fit_results/$1
  [ "$m" = 1 ] && lo=100
  rm -rf $P; mkdir -p $P/HZAmassInde_fTest
  ( cd $F/Background && $S/bin/fTest_ALP_turnOn -i $MV/mA_M$m/ws/run3.root --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root \
      -D $P/HZAmassInde_fTest --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, \
      --mhLow $lo --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 ) > $P/ftest.stdout.log 2>&1
  echo "rc=$?" > $P/.ftest_rc
}
one 1 & one 2 & wait
cat $D/fit_results/*/.ftest_rc > $D/logs/ftest.DONE
