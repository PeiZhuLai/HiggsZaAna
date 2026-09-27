#!/usr/bin/env bash
# Stage 3: once the hadd stage passes GATE 4, regenerate the sideband reweight JSON
# from the FSR-fix ROOT files and compare it against the one the ROOT files were
# actually built with.
#
# The circularity: Parque2Root_BDT.py stamps sideband-reweight branches into every ROOT
# file from sideband_run3_iterative.json, but that JSON is derived from Data/run3.root
# and All_Bkg/run3.root, which are p2root products. The FSR-fix ROOT files were built
# with the OLD JSON (2026-06-03, pre-FSR-fix production).
#
# This only matters if the JSON moved by more than its own statistical noise. The new
# JSON is written to a SEPARATE file; the production one is not touched until the
# comparison has been read.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
M=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA
L=$D/logs_fsrfix/stage3_sideband.log

while [ "$(systemctl --user is-active hza-prepare 2>/dev/null)" = "active" ]; do sleep 60; done
{
echo "[$(date '+%F %T')] prepare finished"
if ! grep -q "VERDICT: PASSED" $D/logs_fsrfix/gate4_report.txt 2>/dev/null; then
  echo "GATE 4 did not pass -- stopping. See logs_fsrfix/gate4_report.txt"
  exit 1
fi
echo "GATE 4 passed, regenerating the sideband reweight JSON from the FSR-fix ROOT"
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
cd $M/scripts
python3 1_make_sideband_reweight.py \
  --output-json $M/reweights/sideband_run3_iterative_fsrfix.json || exit 1
echo
echo "=== comparison: production JSON (pre-FSR-fix) vs newly derived (FSR-fix) ==="
python3 $D/compare_sideband_json.py \
  $M/reweights/sideband_run3_iterative.json \
  $M/reweights/sideband_run3_iterative_fsrfix.json
echo "[$(date '+%F %T')] stage 3 done"
} > $L 2>&1
