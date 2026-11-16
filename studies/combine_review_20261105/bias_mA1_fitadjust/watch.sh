#!/bin/bash
# Watch the mA1 side-study clusters (x0 extension 13463331, fit variants 13463330/2). Exits when no job of them is left in the queue (ALL_DONE), when a job is held
# (HELD, to be handled), or after the hard cap of MAXMIN minutes (TIMEOUT, re-arm). Only condor_q is polled.
W=$(cd "$(dirname "$0")" && pwd); CL="${CL:-13463330 13463331 13463332}"; MAXMIN=${MAXMIN:-240}
for ((i=0; i<MAXMIN; i+=5)); do
  line=""; tot=0; held=0
  for c in $CL; do
    read -r n h r < <(condor_q $c -af JobStatus 2>/dev/null | awk '{n++; if($1==5)h++; if($1==2)r++} END{print n+0, h+0, r+0}')
    line="$line $c:queued=$n,run=$r,held=$h"; tot=$((tot+n)); held=$((held+h))
  done
  echo "$(date '+%F %T')$line" >> $W/logs/watch.txt
  [ $held -gt 0 ] && { echo HELD; exit 0; }
  [ $tot -eq 0 ] && { date '+%F %T' > $W/logs/condor_ALL_DONE; echo ALL_DONE; exit 0; }
  sleep 300
done
echo TIMEOUT
