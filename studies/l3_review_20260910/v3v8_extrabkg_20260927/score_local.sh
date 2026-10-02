#!/usr/bin/env bash
# Local scoring (Parque2Root_BDT.py, same converter + models as run3_bdt_scored_fsrfix):
#   bash score_local.sh <input merged_nominal.parquet> <output .root> [extra converter args]
# Bkg joblist lines of the FSR-fix scoring use ARG_SPLIT=0 -> no --split.
set -uo pipefail
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana || { echo "FATAL: conda activate higgs-alp-ana failed"; exit 1; }
set -u
IN=$1; OUT=$2; shift 2
mkdir -p "$(dirname "$OUT")"
python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py -i "$IN" -o "$OUT" "$@"
