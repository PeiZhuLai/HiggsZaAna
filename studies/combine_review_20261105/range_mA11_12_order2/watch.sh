#!/bin/bash
# Waits for the order-2 bias screen (condor 13433737, 140 chunks: out/o2/*/chunk*.done) and the sideband GOF toys
# (gof/logs/gof_toys.DONE). Completion from markers only. Every 30 min; exits when all done, on a held job, or after 24 h.
# Output: logs/watch.txt, logs/watch.DONE
N=$(cd "$(dirname "$0")" && pwd); L=$N/logs; CL=13433737
rm -f $L/watch.DONE
for i in $(seq 1 48); do
  k=$(find $N/bias/out/o2 -name 'chunk*.done' 2>/dev/null | wc -l)
  g=$([ -f $N/gof/logs/gof_toys.DONE ] && echo done || echo run)
  q=$(condor_q $CL -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q $CL -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') bias $k/140 gof $g | queue[$q] held $held" >> $L/watch.txt
  [ "$k" -ge 140 ] && [ "$g" = done ] && { echo ALL_DONE >> $L/watch.txt; date > $L/watch.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $L/watch.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $L/watch.txt; exit 4
