#!/usr/bin/env bash
# Stage 2 of the FSR-fix chain: hadd the per-era ROOT files into the run3.root files
# the BDT and the fits consume. Runs 2_prepare_rootfile.sh, then GATE 4 (hadd closure).
#
# 2_prepare_rootfile.sh now reads $BASE, which defaults to run3_bdt_inputs_fsrfix.
# Its exit code is not a verdict -- hadd exits 0 while silently dropping trees whose
# branches do not match the FIRST input. GATE 4 is what decides.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
echo "[env] hadd=$(which hadd)  python=$(which python)"
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
bash 2_prepare_rootfile.sh
echo "=== 2_prepare_rootfile.sh raw exit: $? (not a verdict) ==="
echo
python $D/gate4_hadd_closure.py
echo "=== GATE 4 exit: $? ==="
