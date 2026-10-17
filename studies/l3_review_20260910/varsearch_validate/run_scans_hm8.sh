#!/usr/bin/env bash
# Closure + data-GOF scans of the new high-mass model at the 8 anchors not covered in step 3 (mA4,5,6,7,9,10,20,25;
# user request: every point needs a closure). Same recipe as run_scans.sh: DROP=1, cuts of eval/closure_cuts_new_hm8.txt
# (data-sideband count x 2, 1.5, 1, 0.7, 0.5), then the production model at its working point. One lane, NPAR=3.
# Marker: logs/scans_hm8.DONE
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
cd $V; rm -f logs/scans_hm8.DONE
cuts() { awk -v m=$1 '$1==m {$1=""; print}' eval/closure_cuts_new_hm8.txt; }
for m in 4 5 6 7 9 10 20 25; do DROP=1 NPAR=3 bash scan_closure.sh new $m "$(cuts $m)" > logs/scan_new_drop_mA$m.log 2>&1; done
wpof() { python3 -c "import json; print({int(e['mA']): e['MVAcut'] for e in json.load(open('/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json'))['results']}[$1])"; }
pids=()
for m in 4 5 6 7 9 10 20 25; do
  DROP=1 NPAR=1 bash scan_closure.sh prod $m "$(wpof $m)" > logs/scan_prod_drop_mA$m.log 2>&1 & pids+=($!)
  [ ${#pids[@]} -ge 3 ] && { wait ${pids[0]}; pids=("${pids[@]:1}"); }
done
wait
for m in 4 5 6 7 9 10 20 25; do for w in new prod; do echo "== ${w}_drop_mA$m"; cat /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan/${w}_drop_mA$m/summary.txt; done; done > logs/scans_hm8.DONE
