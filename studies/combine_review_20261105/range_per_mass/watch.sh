#!/bin/bash
# Waits for clusters 13409721 (wide-truth, 260 chunks) and 13409722 (highstat consistent ranges, 270 chunks).
# Completion is read from the outputs only: chunk<c>.done for the wide-truth chunks, and for highstat the number of
# BiasFits/*_fits.root per chunk equal to 2 x (truth functions of that mass, read from the moved previous run).
# Checks every 30 min; exits on all done, on any held job, or after 24 h. Status line appended to logs/watch_status.txt.
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/range_per_mass
BD=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal
MASSES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 26 29 30"
for i in $(seq 1 48); do
  a=$(find $R/out -name 'chunk*.done' 2>/dev/null | wc -l)
  b=0
  for m in $MASSES; do
    nf=$(ls $BD/bias_outputs_highstat_variants/mA_${m}_widefit_20261112/chunk0/BiasFits/*_fits.root 2>/dev/null | wc -l)
    for c in $(seq 0 9); do
      n=$(ls $BD/bias_outputs_highstat/mA_$m/chunk$c/BiasFits/*_fits.root 2>/dev/null | wc -l)
      [ "$nf" -gt 0 ] && [ "$n" -eq "$nf" ] && b=$((b+1))
    done
  done
  q=$(condor_q 13409721 13409722 -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  held=$(condor_q 13409721 13409722 -af JobStatus 2>/dev/null | grep -c '^5$')
  echo "$(date '+%F %H:%M') wide-truth $a/260  highstat $b/270  queue[$q] held $held" >> $R/logs/watch_status.txt
  if [ "$a" -eq 260 ] && [ "$b" -eq 270 ]; then echo ALL_DONE >> $R/logs/watch_status.txt; date > $R/logs/ALL.DONE; exit 0; fi
  [ "$held" -gt 0 ] && { echo HELD >> $R/logs/watch_status.txt; exit 3; }
  sleep 1800
done
echo TIMEOUT >> $R/logs/watch_status.txt; exit 4
