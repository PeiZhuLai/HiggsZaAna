#!/bin/bash
# Waits for the jobs of one condor cluster of the 125.38 bias check: returns when every job has left the queue or is held,
# or after a hard limit. Completion marker: "[OK] mA=" lines in the job stdout files of that cluster.
cl="$1"; n="$2"; cap="${3:-14400}"; M=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/bias_mh12538
t0=$(date +%s)
while :; do
  ok=$(cat $M/logs/*.${cl}.*.out 2>/dev/null | grep -c '^\[OK\] mA=')
  q=$(condor_q $cl -af JobStatus 2>/dev/null)
  left=$(echo "$q" | grep -c '[125]'); held=$(echo "$q" | grep -c '^5$')
  echo "$(date +%H:%M) ok=$ok/$n inqueue=$left held=$held"
  if [[ $ok -ge $n || ( $left -eq 0 ) || ( $held -gt 0 && $left -eq $held ) ]]; then break; fi
  (( $(date +%s) - t0 > cap )) && { echo TIMEOUT; break; }
  sleep 300
done
