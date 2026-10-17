#!/bin/bash
# Submit 1_grand_merged.sh as ONE condor job (survives terminal-session / ssh restarts).
# Pass stage switches as env, e.g.
#   RUN_P2ROOT=0 RUN_TREE2WS=0 RUN_SIG_MODEL=0 RUN_BKG=0 bash submit_grand_merged.sh
# 4 cores = the chain's MAX_PARALLEL; 12000 MB keeps 3 GB/core so the schedd
# does not silently add cores. No x509userproxy: the chain needs no grid access.
set -e
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/condor
TS=$(date +%Y%m%d_%H%M%S)
LOGD=$C/logs/$TS; mkdir -p "$LOGD"
CHAIN_LOG=$LOGD/chain.log
envs="CHAIN_LOG=$CHAIN_LOG"
for v in RUN_MERGED_DATA_PRODUCTION RUN_P2ROOT RUN_TREE2WS RUN_SIG_MODEL RUN_BKG \
         RUN_DATACARD RUN_COMBINE RUN_PLOT_LIMITS RUN_BLIND MASSES MERGED_ML_VERSION \
         MERGED_PROD; do
    [ -n "${!v+x}" ] && envs="$envs $v='${!v}'"
done
cat > "$LOGD/job.sub" <<SUB
universe       = vanilla
executable     = $C/run_grand_merged.sh
environment    = "$envs"
getenv         = false
output         = $LOGD/job.out
error          = $LOGD/job.err
log            = $LOGD/job.log
request_cpus   = 4
request_memory = 12000
+JobFlavour    = "workday"
+MaxWallTimeMins = 480
JobBatchName   = "HZa_merged_chain_$TS"
queue 1
SUB
condor_submit "$LOGD/job.sub"
echo "chain log : $CHAIN_LOG"
echo "condor log: $LOGD/job.log"
