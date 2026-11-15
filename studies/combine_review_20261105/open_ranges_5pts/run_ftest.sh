#!/usr/bin/env bash
# Staging fTest with the opened ranges, same command line as shellScripts/bkg/fit_bkg.sh (mhLow 95, GOF toys 100).
# mA19 is the regression point: its ranges are unchanged, so its multipdf must match production.
# Input : /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M<m>/ws/run3.root, Background/bin/fTest_ALP_turnOn (staging build)
# Output: fit_results/<m>/{CMS-HGG_mva_13p6TeV_multipdf.root,HZAmassInde_fTest/,ftest.stdout.log,.done}
# Usage : ./run_ftest.sh <m>            (single start, fit_results/)
#         MS=1 ./run_ftest.sh <m>       (multi-start binary bin/fTest_ALP_turnOn.ms, fit_results_ms/)
#         MS=1 EXTRA_ENV_PDFS=Pow3 ./run_ftest.sh <m>   (one more envelope member, fit_results_ms_extra_Pow3/)
set -o pipefail
m=$1; D=$(cd "$(dirname "$0")" && pwd)
BIN=./bin/fTest_ALP_turnOn; RES=fit_results
if [ "${MS:-0}" = 1 ]; then BIN=./bin/fTest_ALP_turnOn.ms; RES=fit_results_ms; fi
# EXTRA_ENV_PDFS=<names> (patch_extraenv.py) goes to fit_results_ms_extra_<names>/; RES_OVERRIDE sets the directory
if [ -n "${EXTRA_ENV_PDFS:-}" ]; then RES=fit_results_ms_extra_${EXTRA_ENV_PDFS//,/_}; fi
RES=${RES_OVERRIDE:-$RES}; BIN=${BIN_OVERRIDE:-$BIN}
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
O=$D/$RES/$m; mkdir -p $O/HZAmassInde_fTest; rm -f $O/.done
cd $D/Background
$BIN -i /eos/home-p/pelai/HZa/root_MVAcut/data/mA_M$m/ws/run3.root --saveMultiPdf $O/CMS-HGG_mva_13p6TeV_multipdf.root \
  -D $O/HZAmassInde_fTest --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 \
  --gtoys 100 > $O/ftest.stdout.log 2>&1
rc=$?; echo "rc=$rc" > $O/.done
