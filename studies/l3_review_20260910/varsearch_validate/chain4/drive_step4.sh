#!/usr/bin/env bash
# Detached driver of step 4: waits for (a) the signal re-scoring condor cluster(s) and (b) the high-mass closure
# scans (logs/scans_hm8.DONE); resubmits signal samples whose 5 outputs are not all present (up to 3 rounds);
# checks the scored copy (1208 files, no .partial); runs choose_wp.py (stops if a mass has no passing cut); then
# chain4.sh. One status line per poll in logs/drive_step4.status; DONE / STOP / DIED lines in
# ../logs/PENDING_RESULTS.txt (numbers only).
# Usage: drive_step4.sh <rescore cluster id>
CH=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/chain4
V=$(dirname $CH); C=$V/condor; L=$CH/logs; mkdir -p $L
REG=$V/logs/PENDING_RESULTS.txt; ST=$L/drive_step4.status
SBE=/eos/home-p/pelai/HZa/sandbox_varsearch_v1; OUT=$SBE/run3_bdt_scored
ALP_PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
CLUSTERS="$1"
reg() { echo "$(date '+%F %H:%M') $*" >> $REG; }
nsig() { local n=0; for s in $(cat $C/joblist_sig.txt); do for e in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do [ -s $OUT/$s/$e.root ] && n=$((n+1)); done; done; echo $n; }
inq() { local n=0; for c in $CLUSTERS; do n=$((n + $(condor_q $c -af ProcId 2>/dev/null | wc -l))); done; echo $n; }
round=0
while true; do
  q=$(inq); ns=$(nsig); held=0
  for c in $CLUSTERS; do held=$((held + $(condor_q $c -af JobStatus 2>/dev/null | grep -c '^5$'))); done
  sc=$(cat $V/logs/watchdog_hm8.status 2>/dev/null)
  echo "poll $(date '+%F %H:%M:%S') rescore $ns/1190 in queue $q held $held round $round | scans: $sc" > $ST
  if [ "$q" -eq 0 ]; then
    if [ "$ns" -lt 1190 ]; then
      [ $round -ge 3 ] && { reg "STOP  step4 rescore  $ns/1190 after 3 resubmission rounds -> $C/logs"; exit 1; }
      round=$((round+1))
      for s in $(cat $C/joblist_sig.txt); do k=0; for e in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do [ -s $OUT/$s/$e.root ] && k=$((k+1)); done; [ $k -lt 5 ] && echo $s; done > $C/joblist_retry$round.txt
      id=$(cd $C && condor_submit rescore_sig.sub LIST=joblist_retry$round.txt OUTB=$OUT 2>&1 | grep -oP 'submitted to cluster \K[0-9]+')
      [ -n "$id" ] || { reg "STOP  step4 rescore  resubmission round $round failed"; exit 1; }
      CLUSTERS="$CLUSTERS $id"; echo "$(date '+%F %T') round $round: $(wc -l < $C/joblist_retry$round.txt) samples -> cluster $id" >> $L/drive_step4.log
      sleep 300; continue
    fi
    [ -s $V/logs/scans_hm8.DONE ] && break
    grep -q "DIED  hm8\|TIMEOUT hm8" $REG && { reg "STOP  step4  closure scans died (see watchdog_hm8.status)"; exit 1; }
  fi
  sleep 600
done
np=$(find $OUT -name '*.partial' | wc -l); nt=$(find $OUT -name '*.root' -newer $V/models/model_high.pkl | wc -l)
[ "$np" -eq 0 ] && [ "$nt" -eq 1208 ] || { reg "STOP  step4  scored copy: $nt/1208 files, $np partial"; exit 1; }
grep -h "^\[FAILED\]\|GATE FAIL\|READBACK FAIL" $C/logs/*.out > $L/rescore_failures.txt 2>/dev/null
echo "$(date '+%F %T') rescoring complete: 1208 files; failure lines in job logs (retried): $(wc -l < $L/rescore_failures.txt)" >> $L/drive_step4.log
( cd $V && env -i HOME=$HOME PATH=/usr/bin:/bin $ALP_PY choose_wp.py ) > $L/choose_wp.log 2>&1 || { reg "STOP  step4  choose_wp: a mass has no passing cut -> $L/choose_wp.log"; exit 1; }
echo "$(date '+%F %T') working points: $(tr '\n' ';' < $V/eval/wp_chain4.txt)" >> $L/drive_step4.log
bash $CH/chain4.sh > $L/chain4.log 2>&1
if [ -e $L/chain4.DONE ]; then reg "DONE  step4 chain4 (stages 1-10; bias jobs submitted)  -> $L/chain4.log"
else reg "STOP  step4 chain4  $(cat $L/chain4.FAILED 2>/dev/null)  -> $L/chain4.log"; exit 1; fi
