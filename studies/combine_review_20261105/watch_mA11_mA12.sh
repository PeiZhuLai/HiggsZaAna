#!/bin/bash
# Waits for the mA11 deployment checks (condor 13427269 bias 10 chunks, 13427270/13427271 impacts) and the four mA12
# side studies (range_mA12_{floor2,floor1p5,floor3,fixsig2}: driver.DONE and 80 bias chunks each, clusters read from
# bias/logs/condor_clusters.txt). Completion from outputs and markers only. Every 30 min; exits when all done, on a
# held job, on a driver.FAILED, or after 24 h.
# Output: logs/watch_mA11_mA12.txt, logs/watch_mA11_mA12.DONE
ST=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105; L=$ST/logs; mkdir -p $L
C=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine; BD=$C/Checks/Bias_nominal
V="floor2 floor1p5 floor3 fixsig2"
rm -f $L/watch_mA11_mA12.DONE
for i in $(seq 1 48); do
  b=0; im=0
  for c in $(seq 0 9); do n=$(ls $BD/bias_outputs_highstat/mA_11/chunk$c/BiasFits/*_fits.root 2>/dev/null | wc -l); [ "$n" -eq 8 ] && b=$((b+1)); done
  for o in output_impacts_q50_fixidx_prodrange output_impacts_xsec_q50_fixidx_prodrange; do
    [ $C/$o/11_impacts.json -nt $ST/deploy_mA11/logs/stage4.DONE ] && im=$((im+1)); done
  cl="13427269 13427270 13427271"; line="mA11 bias $b/10 impacts $im/2"; alldone=1
  [ "$b" -eq 10 ] && [ "$im" -eq 2 ] || alldone=0
  for v in $V; do
    N=$ST/range_mA12_$v
    [ -f $N/logs/driver.FAILED ] && { echo "$(date '+%F %H:%M') $line | $v FAILED: $(cat $N/logs/driver.FAILED)" >> $L/watch_mA11_mA12.txt; exit 5; }
    k=$(find $N/bias/out -name 'chunk[0-9].done' 2>/dev/null | wc -l)
    d=$([ -f $N/logs/driver.DONE ] && echo done || echo run)
    cl="$cl $(grep -oE 'cluster [0-9]+' $N/bias/logs/condor_clusters.txt 2>/dev/null | awk '{print $2}' | tr '\n' ' ')"
    line="$line | $v $d bias $k/80"
    [ "$d" = done ] && [ "$k" -eq 80 ] || alldone=0
  done
  q=$(condor_q $cl -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q $cl -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') $line | queue[$q] held $held" >> $L/watch_mA11_mA12.txt
  [ "$alldone" -eq 1 ] && { echo ALL_DONE >> $L/watch_mA11_mA12.txt; date > $L/watch_mA11_mA12.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $L/watch_mA11_mA12.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $L/watch_mA11_mA12.txt; exit 4
