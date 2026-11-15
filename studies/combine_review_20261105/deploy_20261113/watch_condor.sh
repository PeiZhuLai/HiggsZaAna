#!/bin/bash
# Waits for clusters 13420172 (bias highstat, 70 chunks), 13420173/13420174 (impacts wilson/xsec, 7 each).
# Completion from the outputs only: per bias chunk, BiasFits/*_fits.root count = 2 x 4 truth functions; per impact,
# <m>_impacts.json newer than logs/stage4.DONE. Every 30 min; exits on all done, any held job, or after 24 h.
# Status appended to logs/watch_status.txt; logs/condor.DONE when complete.
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/deploy_20261113
C=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine; BD=$C/Checks/Bias_nominal
MASSES="9 11 12 13 16 17 20"
for i in $(seq 1 48); do
  b=0; im=0
  for m in $MASSES; do
    for c in $(seq 0 9); do
      n=$(ls $BD/bias_outputs_highstat/mA_$m/chunk$c/BiasFits/*_fits.root 2>/dev/null | wc -l); [ "$n" -eq 8 ] && b=$((b+1))
    done
    for o in output_impacts_q50_fixidx_prodrange output_impacts_xsec_q50_fixidx_prodrange; do
      [ $C/$o/${m}_impacts.json -nt $D/logs/stage4.DONE ] && im=$((im+1))
    done
  done
  q=$(condor_q 13420172 13420173 13420174 -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q 13420172 13420173 13420174 -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') bias $b/70  impacts $im/14  queue[$q] held $held" >> $D/logs/watch_status.txt
  [ "$b" -eq 70 ] && [ "$im" -eq 14 ] && { echo ALL_DONE >> $D/logs/watch_status.txt; date > $D/logs/condor.DONE; exit 0; }
  [ "$held" -gt 0 ] && { echo HELD >> $D/logs/watch_status.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $D/logs/watch_status.txt; exit 4
