#!/bin/bash
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
python3 -c "import sys, numpy, xgboost; print(sys.version.split()[0], numpy.__version__, xgboost.__version__)"
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/bdt_varsearch
mkdir -p $W/test/lcg
OMP_NUM_THREADS=1 NTHREADS=1 python3 $W/search_job.py --samples /eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples --config $W/test/cfg.jsonl --index 0 --out $W/test/lcg
