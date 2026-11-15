#!/usr/bin/env bash
# Highstat bias with the bias fit profiling the workspace (= production limit) parameter ranges at every mass point
# (bias_fit_config.sh default since the 2026-11-12 range consistency fix). Default is a dry run; DRY_RUN=0 does it.
# Masses: all but 21, 27, 28 -- wide-range points without a Bernstein member, where the old bias options set exactly
# the workspace ranges, so their results do not change.
# Old outputs are moved (not deleted) to bias_outputs_highstat_variants/mA_<m>_widefit_20261112.
# Input : Combine/root_t2w (-> root_t2w_xsec)/<m>_Datacard_leptons.root
# Output: Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/, logs Combine/Checks/Bias_nominal/condor_logs_highstat/
set -uo pipefail
R=$(cd "$(dirname "$0")" && pwd); DRY=${DRY_RUN:-1}
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit; BD=$FF/Combine/Checks/Bias_nominal
MASSES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 26 29 30"
ENVI="env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=${KRB5CCNAME:-} XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-}"
do_() { echo "+ $*"; [ "$DRY" = 0 ] && { eval "$@" || { echo "STOP: $*"; exit 1; }; }; true; }
opts=$(bash -c "source $BD/bias_fit_config.sh; unset BIAS_FIT_COMBINE_OPTIONS BIAS_FIT_EXTRA_OPTIONS; resolve_bias_fit_combine_options")
[[ "$opts" == *"--setParameterRanges r=-5,5" ]] || { echo "STOP: bias_fit_config default is not the production-range one"; exit 1; }
for m in $MASSES; do
  src=$BD/bias_outputs_highstat/mA_$m; dst=$BD/bias_outputs_highstat_variants/mA_${m}_widefit_20261112
  [ -e $dst ] && { echo "STOP: $dst exists"; exit 1; }
  [ -d $src ] && do_ mv $src $dst
done
do_ "( cd $FF/shellScripts/bias/Condor && $ENVI bash subjob_bias_highstat.sh $MASSES ) 2>&1 | tee -a $R/logs/condor_submit_highstat.log"
echo "done (DRY_RUN=$DRY)"
