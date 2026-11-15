#!/bin/bash
# Waits for the no-floor bias screen of mA11 (condor 13455274, 40 chunks: bias/out/nf/*/chunk*.done) and the sideband
# GOF toys (gof/logs/gof_toys.DONE). Completion from markers only. Every 30 min; exits when all done, on a held job,
# or after 24 h. Output: logs/watch.txt, logs/watch.DONE
N=$(cd "$(dirname "$0")" && pwd); L=$N/logs; CL=13455274
rm -f $L/watch.DONE
for i in $(seq 1 48); do
  k=$(find $N/bias/out/nf -name 'chunk*.done' 2>/dev/null | wc -l)
  g=$([ -f $N/gof/logs/gof_toys.DONE ] && echo done || echo run)
  q=$(condor_q $CL -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q $CL -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') bias $k/40 gof $g | queue[$q] held $held" >> $L/watch.txt
  [ "$k" -ge 40 ] && [ "$g" = done ] && { echo ALL_DONE >> $L/watch.txt; date > $L/watch.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $L/watch.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $L/watch.txt; exit 4
