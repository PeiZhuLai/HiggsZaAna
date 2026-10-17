#!/usr/bin/env bash
# Re-score every file in joblist_rescore.txt (4 in parallel, lxplus local). Each job writes
# logs/rescore_<sample>_<era>.log; files whose output already exists are skipped. Completion marker:
# logs/rescore.DONE with the number of failures.
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain
O=/eos/home-p/pelai/HZa/sandbox_lowmass_mfree/run3_bdt_scored
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
cd $C; rm -f logs/rescore.DONE
one() { s=$1; e=$2; [ -s $O/$s/$e.root ] && { echo "SKIP $s $e"; return 0; }
  for t in 1 2 3; do $PY rescore_lowmass.py --sample $s --era $e --out-base $O > logs/rescore_${s}_${e}.log 2>&1 && return 0
    grep -q "GATE FAIL\|READBACK FAIL" logs/rescore_${s}_${e}.log && return 1; sleep 20; done; return 1; }
export -f one; export C O PY
cat joblist_rescore.txt | xargs -P4 -L1 bash -c 'one $0 $1 || echo "FAILED $0 $1"' > logs/rescore_driver.log 2>&1
echo "failures=$(grep -c FAILED logs/rescore_driver.log)" > logs/rescore.DONE
