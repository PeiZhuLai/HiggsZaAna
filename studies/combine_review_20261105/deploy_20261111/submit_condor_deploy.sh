#!/usr/bin/env bash
# Condor checks after the local deployment (deploy_plan section 5). Default is a dry run: prints what it would move
# and submit; DRY_RUN=0 does it.
#  bias   : highstat bias for the 11 redeployed masses (10, 15, 21-29), 10 chunks x 1000 toys = 110 jobs,
#           production subjob_bias_highstat.sh (1 core, 3000 MB, tomorrow), workspaces Combine/root_t2w (-> root_t2w_xsec,
#           rebuilt in stage 5). Old outputs moved to bias_outputs_highstat_variants/mA_<m>_before_deploy_20261111.
#  impacts: fixed-function expected impacts, q50 injection, production ranges, IMPACT_SKIP_T2W=1, both cards x 30
#           masses = 60 jobs, 1 core / 3000 MB / IMPACT_PARALLEL=1 (the settings of the 11-09 review run, not the
#           4-core production submit). Old outputs moved to <dir>_bak_deploy_20261111.
# Input : Combine/root_t2w_{xsec,wilson}/<m>_Datacard_leptons.root
# Output: Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/  (logs: Combine/Checks/Bias_nominal/condor_logs_highstat/)
#         Combine/output_impacts_q50_fixidx_prodrange/, Combine/output_impacts_xsec_q50_fixidx_prodrange/
#         (logs: deploy_20261111/logs/impact/), clusters in deploy_20261111/logs/condor_clusters.txt
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd); L=$D/logs; DRY=${DRY_RUN:-1}
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit; C=$FF/Combine
BD=$C/Checks/Bias_nominal; SFX=bak_deploy_20261111
MASSES="10 15 21 22 23 24 25 26 27 28 29"
ENVI="env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=${KRB5CCNAME:-} XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-}"
do_() { echo "+ $*"; [ "$DRY" = 0 ] && { eval "$@" || { echo "STOP: $*"; exit 1; }; }; true; }
[ -f $L/chain.DONE ] || { echo "STOP: local chain not done"; exit 1; }
for m in $MASSES $(seq 1 30); do for k in xsec wilson; do
  [ -s $C/root_t2w_$k/${m}_Datacard_leptons.root ] || { echo "STOP: no workspace $k mA$m"; exit 1; }; done; done

echo "== bias: 110 jobs"
for m in $MASSES; do
  src=$BD/bias_outputs_highstat/mA_$m; dst=$BD/bias_outputs_highstat_variants/mA_${m}_before_deploy_20261111
  [ -e $dst ] && { echo "STOP: $dst exists"; exit 1; }
  [ -d $src ] && do_ mv $src $dst
done
do_ "( cd $FF/shellScripts/bias/Condor && $ENVI bash subjob_bias_highstat.sh $MASSES ) 2>&1 | tee -a $L/condor_submit_bias.log"

echo "== impacts: 60 jobs"
mkdir -p $L/impact
for k in wilson xsec; do
  [ $k = xsec ] && o=$C/output_impacts_xsec_q50_fixidx_prodrange || o=$C/output_impacts_q50_fixidx_prodrange
  [ -e ${o}_$SFX ] && { echo "STOP: ${o}_$SFX exists"; exit 1; }
  [ -d $o ] && do_ mv $o ${o}_$SFX
  [ $k = xsec ] && [ -d ${o}_s1 ] && { [ -e ${o}_s1_$SFX ] && { echo "STOP: ${o}_s1_$SFX exists"; exit 1; }; do_ mv ${o}_s1 ${o}_s1_$SFX; }
  sed -e "s#/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/logs/impact_fixidx/#$L/impact/#" \
      -e "s#batch_name = HZa_impact_fixidx_#batch_name = HZa_impact_deploy_#" \
      $D/../impact_fixidx_$k.submit > $D/impact_deploy_$k.submit
  grep -q "request_cpus = 1" $D/impact_deploy_$k.submit && grep -q "request_memory = 3000" $D/impact_deploy_$k.submit \
    && grep -q "IMPACT_SKIP_T2W=1" $D/impact_deploy_$k.submit && grep -q "$L/impact/" $D/impact_deploy_$k.submit \
    || { echo "STOP: submit file $k not as expected"; exit 1; }
  do_ "$ENVI condor_submit $D/impact_deploy_$k.submit 2>&1 | tee -a $L/condor_submit_impact.log"
done
[ "$DRY" = 0 ] && grep -h 'submitted to cluster' $L/condor_submit_bias.log $L/condor_submit_impact.log | sed "s/^/$(date +%F_%T) /" >> $L/condor_clusters.txt
echo "done (DRY_RUN=$DRY)"
