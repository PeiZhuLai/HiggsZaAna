#!/usr/bin/env bash
# Detached: waits for the m_a 1-3 signal re-scoring with low_swapPt_relE2 (condor cluster $1), then runs chainlow.sh
# (stages 1-10) under a watchdog that stops it if any log exceeds 200 MB. DONE / STOP line in the varsearch registry.
CL=$1
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chainlow_relE2w123
REG=$(dirname $C)/logs/PENDING_RESULTS.txt; ST=$C/logs/drive.status
while [ "$(condor_q $CL -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do echo "$(date '+%F %T') rescore cluster $CL in queue $(condor_q $CL -af ProcId 2>/dev/null | wc -l)" > $ST; sleep 300; done
echo "$(date '+%F %T') chainlow running" > $ST
bash $C/chainlow.sh > $C/logs/chainlow.log 2>&1 &
pid=$!
while kill -0 $pid 2>/dev/null; do
  big=$(find $C/logs $C/flashggFinalFit -name '*.log' -size +200M 2>/dev/null | head -1)
  if [ -n "$big" ]; then
    for p in $(ps -u "$USER" -o pid=,args= | awk '/chainlow_relE2w123\/chainlow\.sh/ {print $1}'); do kill $p; done
    echo "$(date '+%F %H:%M') STOP  chainlow_relE2w123 killed: runaway log $big" >> $REG; exit 1
  fi
  sleep 60
done
if [ -e $C/logs/chain4.DONE ]; then echo "$(date '+%F %H:%M') DONE  chainlow_relE2w123 (mA1-3 full chain, low_relE2_w123 x0.5; bias submitted)  -> $C/logs/chainlow.log" >> $REG
else echo "$(date '+%F %H:%M') STOP  chainlow_relE2w123  $(cat $C/logs/chain4.FAILED 2>/dev/null)" >> $REG; fi
echo "$(date '+%F %T') finished" > $ST
