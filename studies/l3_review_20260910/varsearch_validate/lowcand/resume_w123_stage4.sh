#!/usr/bin/env bash
# Resume the chainlow_relE2w123 chain at stage 4 (Tree2WS) after the data MVA-cut files were redone and checked
# (chainlow_relE2w123/redo_data_mvacut.sh; signal files checked by check_sig_mvacut.sh). Same steps 3-4 as
# followup_w123_full.sh: chain under the drive_chainlow watchdog, then the bias watch.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
REG=$V/logs/PENDING_RESULTS.txt; CH=$V/chainlow_relE2w123; ST=$CH/logs/followup.status
stat() { echo "$(date '+%F %T') $*" > $ST; }
stop() { stat "STOP $*"; echo "$(date '+%F %H:%M') STOP  w123 full chain (resume at 4)  $*" >> $REG; exit 1; }
mv $CH/logs/chainlow.log $CH/logs/chainlow.log.try1_stage3fuse 2>/dev/null
mv $CH/logs/chain4.FAILED $CH/logs/chain4.FAILED.try1_stage3fuse 2>/dev/null
echo "$(date '+%F %H:%M') START w123 full chain resumed at stage 4 -> $ST" >> $REG
stat "chainlow running (resumed at stage 4)"
export START_STAGE=4
bash $CH/drive_chainlow.sh 12997628
[ -e $CH/logs/chain4.DONE ] || stop "chainlow did not finish ($(cat $CH/logs/chain4.FAILED 2>/dev/null))"
BC=$(grep -oP 'submitted to cluster \K[0-9]+' $CH/logs/bias_submit.log | head -1)
[ -n "$BC" ] || stop "no bias cluster in bias_submit.log"
stat "bias cluster $BC"
bash $CH/watch_bias.sh $BC
stat "DONE"
