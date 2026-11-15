#!/bin/bash
# Watchdog of the q50 expected-impact jobs (clusters in impact_q50_clusters.txt): every 30 min writes one
# status line to watch_impact_q50_status.txt; when no job is left in the queue appends DONE (or DIED with
# the missing masses) to PENDING_RESULTS.txt. Numbers only.
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105
C=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine
# clusters are re-read every poll, so resubmissions appended to impact_q50_clusters.txt are tracked
for i in $(seq 1 96); do
  CL=$(grep -oE "cluster [0-9]+|^[0-9]+" $D/impact_q50_clusters.txt | grep -oE "[0-9]+" | tr "\n" " ")
  q=$(condor_q $CL -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  nq=$(condor_q $CL -af JobStatus 2>/dev/null | wc -l)
  nw=$(ls $C/output_impacts_q50/*_impacts.json 2>/dev/null | wc -l)
  nx=$(ls $C/output_impacts_xsec_q50/*_impacts.json 2>/dev/null | wc -l)
  echo "$(date '+%F %H:%M') queue=$nq [$q] json wilson=$nw/30 xsec=$nx/30" >> $D/watch_impact_q50_status.txt
  if [ "$nq" -eq 0 ]; then
    if [ "$nw" -eq 30 ] && [ "$nx" -eq 30 ]; then
      echo "$(date '+%F %H:%M') DONE  impact_q50  -> $C/output_impacts{,_xsec}_q50" >> $D/PENDING_RESULTS.txt
    else
      echo "$(date '+%F %H:%M') DIED  impact_q50  wilson=$nw/30 xsec=$nx/30" >> $D/PENDING_RESULTS.txt
    fi
    exit 0
  fi
  sleep 1800
done
echo "$(date '+%F %H:%M') TIMEOUT impact_q50 after 48h" >> $D/PENDING_RESULTS.txt
