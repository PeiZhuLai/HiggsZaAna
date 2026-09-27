#!/usr/bin/env bash
# Heartbeat monitor for the hza-p2root systemd unit.
# Runs independently of any interactive session (AFS close-to-open means an open log
# is not visible from another node, so we rewrite a small heartbeat file each tick).
# Writes: logs_fsrfix/p2root_heartbeat.txt
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix
OUT=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix
LOG=$D/p2root.log
HB=$D/p2root_heartbeat.txt
EXPECTED=1209
while true; do
  act=$(systemctl --user is-active hza-p2root 2>/dev/null)
  n=$(find "$OUT" -name '*.root' 2>/dev/null | wc -l)
  fail=$( { grep -c "Maximum retries reached" "$LOG" || true; } 2>/dev/null )
  att=$( { grep -c "^Attempt 1:" "$LOG" || true; } 2>/dev/null )
  procs=$( { pgrep -fc "Parque2Root_BDT.py" || true; } 2>/dev/null )
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) unit=$act"
    echo "  root files : $n / $EXPECTED"
    echo "  launched   : $att"
    echo "  hard fails : $fail"
    echo "  running    : $procs python"
    echo "  last log   : $(tail -1 "$LOG" 2>/dev/null | cut -c1-140)"
  } > "$HB"
  [ "$act" != "active" ] && { echo "  UNIT NO LONGER ACTIVE -- monitor exiting" >> "$HB"; break; }
  sleep 600
done
