#!/usr/bin/env bash
# Hand-over: when the 273 files the m_a 1-3 chain needs are patched, stop the 6-process patch run (process group of
# run_patch.sh all), restart the remaining 935 files with 2 processes, and start the production chain (<= 4
# processes), so that at most 6 local processes run at any time. Log: logs/orchestrate.log
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_relE2_20261016
PG=$1
while true; do
  n=0; while read -r s e; do [ -e $D/state/${s}__${e}.ok ] && n=$((n + 1)); done < $D/logs/patch_first.txt
  [ $n -eq 273 ] && break
  [ -e $D/logs/patch.FAILED ] && { echo "$(date '+%F %T') patch FAILED, stop"; exit 1; }
  sleep 60
done
echo "$(date '+%F %T') 273/273 patched; stopping process group $PG"
kill -- -$PG 2>/dev/null; sleep 20
ps -o pid= -g $PG 2>/dev/null | grep -q . && { kill -9 -- -$PG 2>/dev/null; sleep 5; }
rm -rf $D/work/*
setsid bash -c "PATCH_P=2 bash $D/run_patch.sh rest > $D/logs/run_patch_rest.log 2>&1" < /dev/null &
echo "$(date '+%F %T') restarted the remaining files with PATCH_P=2"
setsid bash -c "bash $D/chain_prod_m123.sh > $D/logs/chain_prod_m123.log 2>&1" < /dev/null &
echo "$(date '+%F %T') production chain started"
