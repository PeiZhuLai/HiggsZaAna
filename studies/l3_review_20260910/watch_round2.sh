#!/usr/bin/env bash
# Heartbeat for the round-2 chain. Survives logout; rewrites a small file each tick so
# it is readable from any node (AFS close-to-open makes an open log invisible elsewhere).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
OUT=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix
while true; do
  act=$(systemctl --user is-active hza-round2 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) unit=$act"
    echo "  stage      : $(grep -o '\[[0-9]/6\][^=]*' $D/logs_fsrfix/round2.log 2>/dev/null | tail -1)"
    echo "  root files : $(find $OUT -name '*.root' ! -name 'run3.root' 2>/dev/null | wc -l) (1207 conversions + 6 prepare copies once hadd has run)"
    echo "  launched   : $( { grep -c '^Attempt 1:' $D/logs_fsrfix/round2.log || true; } 2>/dev/null )"
    echo "  hard fails : $( { grep -c 'Maximum retries reached' $D/logs_fsrfix/round2.log || true; } 2>/dev/null )"
    echo "  running    : $( { pgrep -fc 'Parque2Root_BDT.py' || true; } 2>/dev/null ) python"
    echo "  last line  : $(tail -1 $D/logs_fsrfix/round2.log 2>/dev/null | cut -c1-130)"
  } > $D/logs_fsrfix/round2_heartbeat.txt
  [ "$act" != "active" ] && { echo "  UNIT FINISHED -- monitor exiting" >> $D/logs_fsrfix/round2_heartbeat.txt; break; }
  sleep 600
done
