#!/usr/bin/env bash
# Re-score every file of a job list (rescore_selected.py; 4 in parallel, lxplus local). Files whose output exists
# are skipped; transient failures retried twice, gate/read-back failures not. Marker: logs/rescore_<list>.DONE.
# Usage: [OUT_BASE=<dir>] [LOW_MODEL=<NAME of models/model_NAME.pkl, default low>] run_rescore.sh <joblist>
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
O=${OUT_BASE:-/eos/home-p/pelai/HZa/sandbox_varsearch_v1/run3_bdt_scored}; LOWM=${LOW_MODEL:-low}
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
L=$(basename $1 .txt); [ "${LOW_MODEL:-low}" != low ] && L=${L}_${LOW_MODEL}
cd $V; rm -f logs/rescore_$L.DONE
one() { s=$1; e=$2; [ -s $O/$s/$e.root ] && { echo "SKIP $s $e"; return 0; }
  for t in 1 2 3; do env -i HOME=$HOME PATH=/usr/bin:/bin $PY rescore_selected.py --sample $s --era $e --out-base $O --low-model $LOWM > logs/rescore_${LOWM}_${s}_${e}.log 2>&1 && return 0
    grep -q "GATE FAIL\|READBACK FAIL" logs/rescore_${LOWM}_${s}_${e}.log && return 1; sleep 20; done; return 1; }
export -f one; export V O PY LOWM
cat $1 | xargs -P4 -L1 bash -c 'one $0 $1 || echo "FAILED $0 $1"' > logs/rescore_${L}_driver.log 2>&1
echo "failures=$(grep -c FAILED logs/rescore_${L}_driver.log) $(date '+%F %T')" > logs/rescore_$L.DONE
