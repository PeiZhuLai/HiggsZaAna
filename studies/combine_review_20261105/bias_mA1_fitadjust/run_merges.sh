#!/bin/bash
# Merge the mA1 side-study outputs (different from production): x0 extended to 10 chunks, selffit, strat1.
L=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/bias_mA1_fitadjust/logs/merge; mkdir -p $L; rm -f $L/ALL_DONE $L/status.txt
for r in /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal/bias_outputs_inj/x0 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal/bias_outputs_mA1_fitadjust/selffit /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal/bias_outputs_mA1_fitadjust/strat1; do
  t=$(basename $r)
  BASE_DIR=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit INJ_E=0.0 BIAS_HS_ROOT=$r bash /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/bias_mA1_fitadjust/../an_round2/bias_inj/merge_bias_inj.sh 1 > $L/${t}.log 2>&1 &
done
wait
touch $L/ALL_DONE
