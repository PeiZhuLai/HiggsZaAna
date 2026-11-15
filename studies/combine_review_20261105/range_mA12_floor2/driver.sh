#!/usr/bin/env bash
# mA12 2 GeV floor side study, end to end: waits for the production build of mA11 (deploy_mA11 stage 1, at most 2 h),
# then staging build, fTest, envelope parameters, sideband GOF toys (300 per member), bias card + truths, and the
# condor bias screen (4 truths x r=0,q50 x 10 chunks x 1000 toys = 80 jobs, 1 core, 3000 MB).
# Output: see the individual scripts; logs/driver.log, logs/driver.DONE or logs/driver.FAILED; condor cluster in bias/logs/condor_clusters.txt
set -uo pipefail
N=$(cd "$(dirname "$0")" && pwd); L=$N/logs; ST=$(dirname "$N")
die() { echo "$(date '+%F %T') $*" | tee -a $L/driver.log > $L/driver.FAILED; exit 1; }
log() { echo "$(date '+%F %T') $*" >> $L/driver.log; }
rm -f $L/driver.DONE $L/driver.FAILED
for i in $(seq 1 24); do [ -f $ST/deploy_mA11/logs/stage1.DONE ] && break; [ -f $ST/deploy_mA11/logs/chain.FAILED ] && die "deploy_mA11 failed"; sleep 300; done
[ -f $ST/deploy_mA11/logs/stage1.DONE ] || die "deploy_mA11 stage 1 not done after 2 h"
log "build"; bash $N/build_staging.sh > $L/build_staging.stdout.log 2>&1 || die "build_staging"
grep -q 'massALP == 12 ? 2.0' $N/Background/src/PdfModelBuilder.cc || die "floor edit missing"
log "fTest"; bash $N/run_ftest.sh > $L/run_ftest.stdout.log 2>&1; grep -q 'rc=0' $L/ftest.DONE || die "fTest"
log "envelope: $(grep -oE 'pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+' $N/fit_results/12/HZAmassInde_fTest/EnvelopeResults.txt | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' ')"
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
python3 $ST/range_mA11_12_28/dump_params.py $N/fit_results/12/CMS-HGG_mva_13p6TeV_multipdf.root > $L/params_12.txt 2>&1 || log "dump_params failed (not fatal)"
log "bias card + truths"
bash $N/bias/make_card.sh 12 fit_results fl2 > $L/make_card.log 2>&1 || die "make_card"
( cd $N/bias && python3 make_truths.py 12 fl2 ) > $L/make_truths.log 2>&1 || die "make_truths"
( cd $N/bias && python3 gen_submit.py fl2 12 ) >> $L/make_truths.log 2>&1 || die "gen_submit"
grep -q 'request_memory = 3000' $N/bias/bias_screen_fl2_mA12.submit && grep -q 'request_cpus   = 1' $N/bias/bias_screen_fl2_mA12.submit || die "submit file"
env -i HOME=$HOME USER=$USER PATH=$PATH KRB5CCNAME=${KRB5CCNAME:-} XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR:-} condor_submit $N/bias/bias_screen_fl2_mA12.submit > $L/condor_submit.log 2>&1 || die "condor_submit"
grep 'submitted to cluster' $L/condor_submit.log | sed "s/^/$(date +%F_%T) /" >> $N/bias/logs/condor_clusters.txt
log "submitted $(cat $N/bias/logs/condor_clusters.txt)"
log "GOF toys"; bash $N/gof/run_gof_toys.sh > $L/gof.stdout.log 2>&1
[ -f $N/gof/logs/gof_toys.DONE ] || die "GOF toys"
date '+%F %T' > $L/driver.DONE
