#!/usr/bin/env bash
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while true; do
  a=$(systemctl --user is-active hza-datavmc 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) unit=$a"
    echo "  stage      : $(grep -o '\[[0-9]/5\][^=]*' $D/logs_fsrfix/datavmc.log 2>/dev/null | tail -1)"
    echo "  last       : $(tail -1 $D/logs_fsrfix/datavmc.log 2>/dev/null | cut -c1-130)"
    echo "  big2       : $( { condor_q -constraint 'JobBatchName=="p2root_MVAScore_big2"' -format "%d\n" ClusterId 2>/dev/null | wc -l; } ) in queue"
    echo "  dataVmc    : $( { condor_q -constraint 'JobBatchName=="dataVmc_fsrfix"' -format "%d\n" ClusterId 2>/dev/null | wc -l; } ) in queue"
    echo "  scored root: $(find /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix -name '*.root' 2>/dev/null | wc -l) / 1207"
  } > $D/logs_fsrfix/datavmc_heartbeat.txt
  [ "$a" != "active" ] && { echo "  FINISHED" >> $D/logs_fsrfix/datavmc_heartbeat.txt; break; }
  sleep 600
done
