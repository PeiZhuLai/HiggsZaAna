#!/usr/bin/env bash
# Watch one bias-screen batch: exit when all chunk .done markers exist (DONE) or any job of the cluster is held
# (HELD, to be handled), hard cap 6 h. Writes one status line every 5 min to watch_<tag>.txt.
# Usage : watch.sh <tag> <cluster> <expected chunks>
B=$(cd "$(dirname "$0")" && pwd); tag=$1; cl=$2; n=$3; t0=$(date +%s); st=$B/watch_${tag}.txt
while true; do
  d=$(find $B/out/$tag -name 'chunk*.done' 2>/dev/null | wc -l)
  q=$(condor_q -constraint "ClusterId==$cl" -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  h=$(condor_q -constraint "ClusterId==$cl && JobStatus==5" -af ProcId 2>/dev/null | wc -l)
  echo "$(date +%F_%T) done=$d/$n queue[$q] held=$h" >> $st
  [ "$d" -ge "$n" ] && { echo "DONE $tag $d/$n"; exit 0; }
  [ "$h" -gt 0 ] && { echo "HELD $tag $h jobs"; exit 2; }
  [ $(( $(date +%s) - t0 )) -gt 21600 ] && { echo "TIMEOUT $tag $d/$n"; exit 1; }
  sleep 300
done
