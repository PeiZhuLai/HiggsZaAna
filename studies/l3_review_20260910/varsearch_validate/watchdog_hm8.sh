#!/usr/bin/env bash
# Detached watchdog for run_scans_hm8.sh: status line per poll in logs/watchdog_hm8.status; DONE / DIED / TIMEOUT
# line in logs/PENDING_RESULTS.txt. Numbers only.
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
REG=$V/logs/PENDING_RESULTS.txt; ST=$V/logs/watchdog_hm8.status; B=/eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan
for _ in $(seq 1 360); do
  nn=0; np=0
  for m in 4 5 6 7 9 10 20 25; do
    nn=$((nn + $(ls $B/new_drop_mA$m/cut_*/result.txt 2>/dev/null | wc -l))); np=$((np + $(ls $B/prod_drop_mA$m/cut_*/result.txt 2>/dev/null | wc -l)))
  done
  al=$(ps -u "$USER" -o cmd= --no-headers | grep -c '[r]un_scans_hm8.sh' || true)
  echo "poll $(date '+%F %H:%M:%S') new cuts $nn/40 | prod $np/8 | alive=$al" > $ST
  if [ -s $V/logs/scans_hm8.DONE ]; then echo "$(date '+%F %H:%M') DONE  hm8 closure scans  -> $V/logs/scans_hm8.DONE" >> $REG; exit 0; fi
  if [ "$al" -eq 0 ]; then echo "$(date '+%F %H:%M') DIED  hm8 closure scans  new $nn/40 prod $np/8" >> $REG; exit 1; fi
  sleep 60
done
echo "$(date '+%F %H:%M') TIMEOUT hm8 watchdog  new $nn/40 prod $np/8" >> $REG
