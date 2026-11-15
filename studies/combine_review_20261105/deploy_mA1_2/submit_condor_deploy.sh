#!/usr/bin/env bash
# Condor checks after the local deployment of mA1, 2 (chain_deploy.sh stages 1-4).
# Default is a dry run: prints what it would move and submit; DRY_RUN=0 does it.
#  bias   : highstat bias, 10 chunks x 1000 toys per mass = 20 jobs, production subjob_bias_highstat.sh (1 core,
#           3000 MB), workspaces Combine/root_t2w (-> root_t2w_xsec, rebuilt in stage 4), bias fit with the production
#           ranges (bias_fit_config.sh default). The outputs of the 11-12 run with the per-mass windows are moved to
#           bias_outputs_highstat_variants/mA_<m>_narrow_20261112.
#  impacts: fixed-function expected impacts, q50 injection, production ranges, both cards x 2 masses = 4 jobs,
#           same submit settings as deploy_20261111 (1 core / 3000 MB). The 2 masses' old files are moved to
#           <dir>_bak_deploy_mA1_2/.
# Input : Combine/root_t2w_{xsec,wilson}/<m>_Datacard_leptons.root
# Output: Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_<m>/  (logs: Combine/Checks/Bias_nominal/condor_logs_highstat/)
#         Combine/output_impacts_q50_fixidx_prodrange/, Combine/output_impacts_xsec_q50_fixidx_prodrange/
#         (logs: deploy_mA1_2/logs/impact/), clusters in deploy_mA1_2/logs/condor_clusters.txt
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd); L=$D/logs; DRY=${DRY_RUN:-1}
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit; C=$FF/Combine
BD=$C/Checks/Bias_nominal; SFX=bak_deploy_mA1_2
MASSES="1 2"
ENVI="env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=${KRB5CCNAME:-} XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-}"
do_() { echo "+ $*"; [ "$DRY" = 0 ] && { eval "$@" || { echo "STOP: $*"; exit 1; }; }; true; }
[ -f $L/chain.DONE ] || { echo "STOP: local chain not done"; exit 1; }
for m in $MASSES; do for k in xsec wilson; do
  [ $C/root_t2w_$k/${m}_Datacard_leptons.root -nt $L/stage3.DONE ] || { echo "STOP: workspace $k mA$m not rebuilt"; exit 1; }; done; done

echo "== bias: 20 jobs"
for m in $MASSES; do
  src=$BD/bias_outputs_highstat/mA_$m; dst=$BD/bias_outputs_highstat_variants/mA_${m}_narrow_20261112
  [ -e $dst ] && { echo "STOP: $dst exists"; exit 1; }
  [ -d $src ] && do_ mv $src $dst
done
do_ "( cd $FF/shellScripts/bias/Condor && $ENVI bash subjob_bias_highstat.sh $MASSES ) 2>&1 | tee -a $L/condor_submit_bias.log"

echo "== impacts: 4 jobs"
mkdir -p $L/impact
for k in wilson xsec; do
  [ $k = xsec ] && o=$C/output_impacts_xsec_q50_fixidx_prodrange || o=$C/output_impacts_q50_fixidx_prodrange
  [ -e ${o}_$SFX ] && { echo "STOP: ${o}_$SFX exists"; exit 1; }
  do_ mkdir -p ${o}_$SFX
  for m in $MASSES; do for f in ${m}_Expected_Impacts.pdf ${m}_impacts.json; do
    [ -e $o/$f ] && do_ mv $o/$f ${o}_$SFX/
  done; done
  sed -e "s#$D/../deploy_20261111/logs/impact/#$L/impact/#" -e "s#/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/deploy_20261111/logs/impact/#$L/impact/#" \
      -e "s#batch_name = HZa_impact_deploy_#batch_name = HZa_impact_deployA12_#" \
      -e "s#^queue mA in (.*)#queue mA in ($MASSES)#" \
      $D/../deploy_20261111/impact_deploy_$k.submit > $D/impact_deploy_$k.submit
  grep -q "request_cpus = 1" $D/impact_deploy_$k.submit && grep -q "request_memory = 3000" $D/impact_deploy_$k.submit \
    && grep -q "IMPACT_SKIP_T2W=1" $D/impact_deploy_$k.submit && grep -q "$L/impact/" $D/impact_deploy_$k.submit \
    && grep -qx "queue mA in ($MASSES)" $D/impact_deploy_$k.submit && ! grep -q deploy_20261111 $D/impact_deploy_$k.submit \
    || { echo "STOP: submit file $k not as expected"; exit 1; }
  do_ "$ENVI condor_submit $D/impact_deploy_$k.submit 2>&1 | tee -a $L/condor_submit_impact.log"
done
[ "$DRY" = 0 ] && grep -h 'submitted to cluster' $L/condor_submit_bias.log $L/condor_submit_impact.log | sed "s/^/$(date +%F_%T) /" >> $L/condor_clusters.txt
echo "done (DRY_RUN=$DRY)"
