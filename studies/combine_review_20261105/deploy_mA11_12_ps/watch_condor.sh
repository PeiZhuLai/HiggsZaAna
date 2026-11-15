#!/bin/bash
# Waits for the condor checks of the mA11/mA12 deployment: highstat bias (cluster 13446820, 20 jobs) and expected
# impacts (13446821 wilson, 13446822 xsec, 2 jobs each). Done when no job of these clusters is left in the queue and
# every expected output exists: 8 fit files per bias chunk (4 truths x 2 splits) and the 4 impact json files newer
# than the submission. Every 30 min; exits on a held job, or after 24 h.
# Output: logs/watch_condor.txt, logs/watch_condor.DONE
D=$(cd "$(dirname "$0")" && pwd); L=$D/logs
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit; C=$FF/Combine
HS=$C/Checks/Bias_nominal/bias_outputs_highstat; CL="13446820 13446821 13446822"
REF=$L/condor_clusters.txt
rm -f $L/watch_condor.DONE
for i in $(seq 1 48); do
  q=$(condor_q $CL -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  n=$(condor_q $CL -af JobStatus 2>/dev/null | wc -l)
  held=$(condor_q $CL -af JobStatus 2>/dev/null | grep -c '^5$')
  nb=0; for m in 11 12; do for c in $(seq 0 9); do
    k=$(ls $HS/mA_$m/chunk$c/BiasFits/biasStudy_*_split*_fits.root 2>/dev/null | wc -l); [ "$k" -ge 8 ] && nb=$((nb+1)); done; done
  ni=0; for o in $C/output_impacts_q50_fixidx_prodrange $C/output_impacts_xsec_q50_fixidx_prodrange; do for m in 11 12; do
    [ $o/${m}_impacts.json -nt $REF ] && ni=$((ni+1)); done; done
  echo "$(date '+%F %H:%M') queue[$q] held $held | bias chunks $nb/20 impacts $ni/4" >> $L/watch_condor.txt
  [ "$n" -eq 0 ] && [ "$nb" -ge 20 ] && [ "$ni" -ge 4 ] && { echo ALL_DONE >> $L/watch_condor.txt; date > $L/watch_condor.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $L/watch_condor.txt; exit 3; }
  [ "$n" -eq 0 ] && { echo "QUEUE_EMPTY_OUTPUT_MISSING" >> $L/watch_condor.txt; exit 5; }
  sleep 1800
done
echo TIMEOUT >> $L/watch_condor.txt; exit 4
