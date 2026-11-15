#!/usr/bin/env bash
# Verify patch_databinned.py with bin/fTest_ALP_turnOn.msd at mA10 (2 fTests in parallel):
#  A MHLOW=100: data_obs = histogram of the 101 fit-window events; every variable as fit_results_ms_lo100 (old binary)
#  B MHLOW=95 : whole workspace identical to fit_results_ms (old binary), the unchanged path
# Output: fit_results_msd_lo100/10/, fit_results_msd/10/, verify_databinned.txt
D=$(cd "$(dirname "$0")" && pwd); M=CMS-HGG_mva_13p6TeV_multipdf.root
MS=1 MHLOW=100 BIN_OVERRIDE=./bin/fTest_ALP_turnOn.msd RES_OVERRIDE=fit_results_msd_lo100 $D/run_ftest_v3.sh 10 &
MS=1 MHLOW=95  BIN_OVERRIDE=./bin/fTest_ALP_turnOn.msd RES_OVERRIDE=fit_results_msd       $D/run_ftest_v3.sh 10 &
wait
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
{ echo "A mhLow 100: $(cat $D/fit_results_msd_lo100/10/.done)"
  python3 $D/compare_ws.py $D/fit_results_ms_lo100/10/$M $D/fit_results_msd_lo100/10/$M --data-may-differ --expect-data-from-unbinned
  echo "B mhLow 95: $(cat $D/fit_results_msd/10/.done)"
  python3 $D/compare_ws.py $D/fit_results_ms/10/$M $D/fit_results_msd/10/$M
} > $D/verify_databinned.txt 2>&1
cat $D/verify_databinned.txt
