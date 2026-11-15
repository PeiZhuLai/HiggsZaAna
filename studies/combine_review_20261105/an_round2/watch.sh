#!/bin/bash
# Watch the round-2 condor clusters. Exits when no job of them is left in the queue (ALL_DONE), when a job is held
# (HELD, to be handled), or after the hard cap of MAXMIN minutes (TIMEOUT, re-arm). Only condor_q is polled.
W=$(cd "$(dirname "$0")" && pwd); CL="13460619 13460620 13460621"; MAXMIN=${MAXMIN:-90}
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
