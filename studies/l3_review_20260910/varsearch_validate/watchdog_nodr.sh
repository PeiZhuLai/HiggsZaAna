#!/usr/bin/env bash
# Detached watchdog for run_nodr.sh: one status line per poll in logs/watchdog_nodr.status, DONE / DIED / TIMEOUT
# line in logs/PENDING_RESULTS.txt. Numbers only.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
REG=$V/logs/PENDING_RESULTS.txt; ST=$V/logs/watchdog_nodr.status; B=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
for _ in $(seq 1 240); do
  ns=$(ls $B/closure_scan/nodr_drop_mA*/cut_*/result.txt 2>/dev/null | wc -l)
  nb=$(ls $B/closure_boot/nodr_drop_mA1_cut*/rep_*/result.txt 2>/dev/null | wc -l)
  al=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_nodr.sh' || true)
  echo "poll $(date '+%F %H:%M:%S') scan cuts $ns/15 | boot replicas $nb/60 | alive=$al" > $ST
  if [ -e $V/logs/nodr.DONE ]; then echo "$(date '+%F %H:%M') DONE  nodr closure+bootstrap  -> $V/logs/nodr.DONE" >> $REG; exit 0; fi
  if [ "$al" -eq 0 ]; then echo "$(date '+%F %H:%M') DIED  nodr closure+bootstrap  scan $ns/15 boot $nb/60" >> $REG; exit 1; fi
  sleep 60
done
echo "$(date '+%F %H:%M') TIMEOUT nodr watchdog  scan $ns/15 boot $nb/60" >> $REG
