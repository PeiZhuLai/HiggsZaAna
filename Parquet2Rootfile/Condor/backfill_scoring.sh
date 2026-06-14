#!/bin/bash
# Held-tolerant scoring backfill (bypasses 1_grand.sh, which aborts on any held).
# Workers write the .root DIRECTLY to EOS; a job held on the AFS .out transfer
# produced NO root, so we just resubmit whatever's still missing on EOS until 0.
# condor_q failures (flaky schedd) are treated as "retry", never as "drained".
set -uo pipefail
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor

BATCH='p2root_MVAScore'
HELD_C="JobBatchName==\"$BATCH\" && JobStatus==5"
ACTIVE_C="JobBatchName==\"$BATCH\" && (JobStatus==1 || JobStatus==2)"

cq() { # robust condor_q -af count of a constraint; echoes -1 on query failure, else exact count
  local out rc
  out=$(condor_q -constraint "$1" -af ClusterId 2>/dev/null); rc=$?
  [ $rc -ne 0 ] && { echo -1; return; }
  [ -z "$out" ] && { echo 0; return; }
  printf '%s\n' "$out" | grep -c .
}

count_missing() {
  python3 1_make_joblist.py >/dev/null 2>&1   # full 1209 -> joblist.tsv
  : > joblist_missing.tsv; local m=0
  while IFS=$'\t' read -r inp out corr split; do
    [ -f "$out" ] || { printf '%s\t%s\t%s\t%s\n' "$inp" "$out" "$corr" "$split" >> joblist_missing.tsv; m=$((m+1)); }
  done < joblist.tsv
  echo "$m"
}

for attempt in $(seq 1 40); do
  echo "[$(date +%H:%M:%S)] === attempt $attempt: drain active jobs ==="
  # Drain: wait until no idle/running batch jobs; release-then-rm held (they made no root).
  while true; do
    h=$(cq "$HELD_C"); [ "$h" -gt 0 ] 2>/dev/null && condor_rm -constraint "$HELD_C" >/dev/null 2>&1
    n=$(cq "$ACTIVE_C")
    if [ "$n" = "-1" ]; then echo "  [warn] condor_q failed; retry in 30s"; sleep 30; continue; fi
    echo "  active=$n held=$h"
    [ "$n" -eq 0 ] && break
    sleep 30
  done

  miss=$(count_missing)
  echo "[$(date +%H:%M:%S)] missing roots = $miss"
  if [ "$miss" -eq 0 ]; then echo "ALL_SCORED"; exit 0; fi

  cp joblist_missing.tsv joblist.tsv
  echo "[$(date +%H:%M:%S)] submitting $miss missing jobs"
  condor_submit 2_submit.sub 2>&1 | tail -2
  sleep 20
done

echo "BACKFILL_INCOMPLETE after 40 attempts; missing=$(count_missing)"
exit 1
