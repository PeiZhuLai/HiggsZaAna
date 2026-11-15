#!/usr/bin/env bash
# Watch the deployment condor clusters (bias 13403435, impacts 13403436/13403437) every 10 minutes, hard cap 24 h.
# Exits 0 when no job of the clusters is left in the queue, 2 when a job is held (for inspection), 3 at the cap.
# Output: logs/watch_condor.txt (one status line per poll)
L=$(cd "$(dirname "$0")" && pwd)/logs; CL="13403435 13403436 13403437"; end=$(( $(date +%s) + 86400 ))
while [ "$(date +%s)" -lt "$end" ]; do
  line="$(date '+%F %T')"; left=0; held=0
  for c in $CL; do
    s=$(condor_q $c -af JobStatus 2>/dev/null | sort | uniq -c | awk '{printf "%s:%s ", $2, $1}')
    n=$(condor_q $c -af JobStatus 2>/dev/null | wc -l); h=$(condor_q $c -af JobStatus 2>/dev/null | grep -c '^5$')
    left=$((left + n)); held=$((held + h)); line="$line | $c ${s:-empty}"
  done
  echo "$line" >> $L/watch_condor.txt
  [ "$held" -gt 0 ] && { echo "HELD $held" >> $L/watch_condor.txt; exit 2; }
  [ "$left" -eq 0 ] && { echo "ALL DONE" >> $L/watch_condor.txt; exit 0; }
  sleep 600
done
echo "CAP reached" >> $L/watch_condor.txt; exit 3
