#!/usr/bin/env bash
# Closure + data-GOF scans (step 3 of the validation of the variable-search BDTs), pseudo-data without the BDT
# training DY events (DROP=1). Two lanes x NPAR 3 = 6 processes (lxplus foreground limit).
#   lane A: new model mA1, 2, 3 at the cuts of eval/closure_cuts_new.txt (data-sideband count x 2, 1.5, 1, 0.7, 0.5)
#   lane B: new model mA8, 15, 30; then the production model at its working point (same DROP) for all six masses
# Marker: logs/scans.DONE
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
cd $V; rm -f logs/scans.DONE
cuts() { awk -v m=$1 '$1==m {$1=""; print}' eval/closure_cuts_new.txt; }
laneA() { for m in 1 2 3; do DROP=1 NPAR=3 bash scan_closure.sh new $m "$(cuts $m)" > logs/scan_new_drop_mA$m.log 2>&1; done; }
laneB() { for m in 8 15 30; do DROP=1 NPAR=3 bash scan_closure.sh new $m "$(cuts $m)" > logs/scan_new_drop_mA$m.log 2>&1; done
  for m in 1 2 3 8 15 30; do
    wp=$(python3 -c "import json; print({int(e['mA']): e['MVAcut'] for e in json.load(open('/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json'))['results']}[$m])")
    DROP=1 NPAR=1 bash scan_closure.sh prod $m "$wp" > logs/scan_prod_drop_mA$m.log 2>&1
  done; }
laneA & laneB & wait
cat /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_scan/*/summary.txt > logs/scans.DONE
