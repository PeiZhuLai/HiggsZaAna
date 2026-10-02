#!/usr/bin/env bash
# p2root with the same converter/env as run_p2root_fsrfix.sh (higgs-alp-ana; hza_ana lacks xgboost
# and breaks ROOT after pandas):  bash p2root_local.sh <merged_nominal.parquet> <out.root> [args]
set -uo pipefail
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana || { echo "FATAL: conda activate higgs-alp-ana failed"; exit 1; }
python -c "
import pandas, uproot
from ROOT import Math, TVector2, TVector3, TLorentzVector
from xgboost import XGBClassifier
" || { echo "FATAL: preflight imports failed"; exit 1; }
set -u
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
IN=$1; OUT=$2; shift 2
mkdir -p "$(dirname "$OUT")"
python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py -i "$IN" -o "$OUT" "$@"
