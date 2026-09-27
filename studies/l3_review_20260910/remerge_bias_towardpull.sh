#!/usr/bin/env bash
# Re-merge the high-stat bias of every mass point except m_a = 20 (being rerun) with the new pull
# definition in RunBiasStudy.py (error facing the truth, 2026-09-26). No new toys: the merge only
# re-hadds the per-chunk fits and re-runs the Gaussian fit + plot_bias.py.
# plots_bias (flat, AN source) is moved aside first so no symmetric-error plot survives.
set -uo pipefail
L=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
BN=$F/Combine/Checks/Bias_nominal
MASSES=$(seq 1 30 | grep -vx 20 | tr '\n' ' ')
ROOT_DATACARD_PATH=$F/Combine/root_t2w_fsrfix_20260925 bash $F/shellScripts/bias/Condor/merge_bias_highstat.sh $MASSES > $L/remerge_towardpull.log 2>&1
echo "merge rc=$?"
[ -d $BN/plots_bias ] && mv $BN/plots_bias $BN/plots_bias_symunc_20260926
BASE_DIR=$F bash $F/shellScripts/bias/Condor/collect_bias_results.sh
n=0; for m in $MASSES; do j=$BN/bias_outputs_highstat/mA_$m/merged/BiasJson/${m}_gaussfit.json; [ -s $j ] && [ $j -nt $BN/RunBiasStudy.py ] && n=$((n+1)); done
echo "fresh gaussfit json (toward-truth pulls): $n/29"
