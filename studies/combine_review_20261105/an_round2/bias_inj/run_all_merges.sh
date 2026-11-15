#!/bin/bash
# Merge the injection bias study (0, 2, 3 x q50) for mA = 1, 10, 20, 30, one (injection, mass) at a time,
# with the injected signal of each point (values from bias_inj.submit). Writes logs/bias_inj_merge/<tag>_mA<m>.log
# and the marker logs/bias_inj_merge/ALL_DONE when finished.
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/an_round2
BD=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal
L=$W/logs/bias_inj_merge
mkdir -p "$L"; rm -f "$L/ALL_DONE" "$L/status.txt"
declare -A E=(
  [x0_1]=0.0 [x0_10]=0.0 [x0_20]=0.0 [x0_30]=0.0
  [x2_1]=0.6094 [x2_10]=0.0383 [x2_20]=0.0836 [x2_30]=0.0891
  [x3_1]=0.9141 [x3_10]=0.0574 [x3_20]=0.1254 [x3_30]=0.1336
)
for tag in x0 x2 x3; do
  for m in 1 10 20 30; do
    BASE_DIR=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit INJ_E=${E[${tag}_${m}]} BIAS_HS_ROOT=$BD/bias_outputs_inj/$tag \
      bash $W/bias_inj/merge_bias_inj.sh $m > "$L/${tag}_mA${m}.log" 2>&1
    echo "$tag mA$m rc=$?" >> "$L/status.txt"
  done
done
touch "$L/ALL_DONE"
