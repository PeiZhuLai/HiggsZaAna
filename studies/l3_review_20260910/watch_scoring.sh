#!/usr/bin/env bash
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while true; do
  a1=$(systemctl --user is-active hza-trainlh 2>/dev/null)
  a2=$(systemctl --user is-active hza-scoring 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) trainlh=$a1 scoring=$a2"
    echo "  train stage : $(grep -o '\[[0-9]/4\][^=]*' $D/logs_fsrfix/train_lowhigh.log 2>/dev/null | tail -1)"
    echo "  train last  : $(tail -1 $D/logs_fsrfix/train_lowhigh.log 2>/dev/null | cut -c1-120)"
    echo "  score stage : $(grep -o '\[[0-9]/5\][^=]*' $D/logs_fsrfix/scoring.log 2>/dev/null | tail -1)"
    echo "  score last  : $(tail -1 $D/logs_fsrfix/scoring.log 2>/dev/null | cut -c1-120)"
    echo "  condor      : $( { condor_q -constraint 'JobBatchName=="p2root_MVAScore" || JobBatchName=="p2root_MVAScore_TEST"' -format "%d\n" ClusterId 2>/dev/null | wc -l; } ) in queue"
    echo "  scored root : $(find /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix -name '*.root' 2>/dev/null | wc -l) / 1207"
  } > $D/logs_fsrfix/scoring_heartbeat.txt
  [ "$a1" != "active" ] && [ "$a2" != "active" ] && { echo "  BOTH FINISHED" >> $D/logs_fsrfix/scoring_heartbeat.txt; break; }
  sleep 600
done
