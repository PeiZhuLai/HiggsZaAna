#!/usr/bin/env bash
# Detached watchdog for the step-3 closure scans (run_scans.sh -> logs/scans.DONE) and the mA1 bootstrap
# (run_boot.sh -> logs/boot.DONE). Writes one line per poll to logs/watchdog.status and DONE / DIED lines to
# logs/PENDING_RESULTS.txt. Numbers only, no interpretation.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
REG=$V/logs/PENDING_RESULTS.txt; ST=$V/logs/watchdog.status; B=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
sdone=0
for _ in $(seq 1 480); do
  ns=$(ls $B/closure_scan/*/summary.txt 2>/dev/null | wc -l)
  nb=$(ls $B/closure_boot/*/rep_*/result.txt 2>/dev/null | wc -l)
  as=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_scans.sh' || true)
  ab=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_boot.sh' || true)
  echo "poll $(date '+%F %H:%M:%S') scan summaries $ns/12 alive=$as | boot replicas $nb/90 alive=$ab" > $ST
  if [ $sdone = 0 ] && [ -e $V/logs/scans.DONE ]; then echo "$(date '+%F %H:%M') DONE  closure scans  -> $V/logs/scans.DONE" >> $REG; sdone=1; fi
  if [ $sdone = 0 ] && [ "$as" -eq 0 ]; then echo "$(date '+%F %H:%M') DIED  closure scans  after $ns/12 summaries" >> $REG; sdone=1; fi
  if [ -e $V/logs/boot.DONE ]; then echo "$(date '+%F %H:%M') DONE  mA1 bootstrap  -> $V/logs/boot.DONE" >> $REG; exit 0; fi
  if [ "$ab" -eq 0 ]; then echo "$(date '+%F %H:%M') DIED  mA1 bootstrap  after $nb/90 replicas" >> $REG; exit 1; fi
  sleep 60
done
echo "$(date '+%F %H:%M') TIMEOUT watchdog  scans $ns/12, boot $nb/90" >> $REG
