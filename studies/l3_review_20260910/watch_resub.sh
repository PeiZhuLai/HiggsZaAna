#!/usr/bin/env bash
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while true; do
  a=$(systemctl --user is-active hza-resub 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) unit=$a"
    echo "  stage      : $(grep -o '\[[0-9]/4\][^=]*' $D/logs_fsrfix/scoring_resub.log 2>/dev/null | tail -1)"
    echo "  last       : $(tail -1 $D/logs_fsrfix/scoring_resub.log 2>/dev/null | cut -c1-130)"
    echo "  condor     : $( { condor_q -constraint 'JobBatchName=="p2root_MVAScore_resub" || JobBatchName=="p2root_MVAScore_TEST"' -format "%d\n" ClusterId 2>/dev/null | wc -l; } ) in queue"
    echo "  scored root: $(find /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix -name '*.root' 2>/dev/null | wc -l) / 1207"
  } > $D/logs_fsrfix/resub_heartbeat.txt
  [ "$a" != "active" ] && { echo "  FINISHED" >> $D/logs_fsrfix/resub_heartbeat.txt; break; }
  sleep 600
done
