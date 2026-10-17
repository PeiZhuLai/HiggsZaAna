#!/bin/bash
# Local check of the new gates A/B on the stored production scores (no training).
set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh; set -u
D=$(dirname $(readlink -f $0)); cd $D/../..
python3 search_job.py --samples /eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples --config $D/config.jsonl --index $1 --out $D/out
