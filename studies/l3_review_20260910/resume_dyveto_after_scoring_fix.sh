#!/usr/bin/env bash
# 2026-09-29 19:5x: resume the dyveto chain after the scoring-retry fix.
#  1. wait for the two big scoring retries (12773003.1194 DYG 2024 at 8 GB, 12778383.0 Data 2024 at 8 GB/workday)
#  2. chain stage D (SKIP_FULL_SCORING=1: only retries + GATE 5) and E
#  3. launch the downstream driver only once the dataVmc code changes are in (flag file)
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
while [ "$(condor_q 12773003.1194 12778383 -format '%d\n' ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  for id in $(condor_q 12773003.1194 12778383 -constraint 'JobStatus==5 && RequestMemory < 16000' -af:j ClusterId 2>/dev/null | sed 's/ .*//'); do
    echo "$(date '+%F %T') $id held -> RequestMemory 16000 + release" >> $L/dyveto_retrain_heartbeat.txt
    condor_qedit $id RequestMemory 16000 >/dev/null; condor_release $id >/dev/null
  done
  sleep 300
done
echo "$(date '+%F %T') scoring retries left the queue; resuming chain at D (SKIP_FULL_SCORING=1)" >> $L/dyveto_retrain_heartbeat.txt
START_STAGE=D SKIP_FULL_SCORING=1 bash $D/chain_dyveto_retrain.sh > $L/dyveto_retrain.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_retrain.log || exit 1
while [ ! -f $L/dyveto_datavmc_code_ready ]; do sleep 120; done
systemd-run --user --unit=hza-dyveto-driver --collect -p WorkingDirectory=$D bash -c "bash $D/chain_dyveto_driver.sh > $L/dyveto_driver.log 2>&1"
