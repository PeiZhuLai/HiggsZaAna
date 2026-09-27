#!/usr/bin/env bash
# Re-arm stage 3 (sideband JSON regeneration + comparison) once the signal reconversion
# and both hadd gates have passed. stage 3 fail-closed on the first attempt because
# GATE 4 had caught the dropped 2024 trees -- correct behavior, it just needs rerunning.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
M=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA
L=$D/logs_fsrfix/stage3_sideband.log
while [ "$(systemctl --user is-active hza-refix 2>/dev/null)" = "active" ]; do sleep 60; done
{
echo "[$(date '+%F %T')] refix finished"
grep -q "VERDICT: PASSED" $D/logs_fsrfix/gate4_report.txt 2>/dev/null || { echo "GATE 4 still failing -- stop"; exit 1; }
grep -q "VERDICT: PASSED" $D/logs_fsrfix/gate4b_report.txt 2>/dev/null || { echo "GATE 4b still failing -- stop"; exit 1; }
echo "both gates passed; regenerating the sideband reweight JSON from the FSR-fix ROOT"
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
cd $M/scripts
python3 1_make_sideband_reweight.py --output-json $M/reweights/sideband_run3_iterative_fsrfix.json || exit 1
echo
echo "=== comparison: production JSON (pre-FSR-fix) vs newly derived (FSR-fix) ==="
python3 $D/compare_sideband_json.py \
  $M/reweights/sideband_run3_iterative.json \
  $M/reweights/sideband_run3_iterative_fsrfix.json
echo "[$(date '+%F %T')] stage 3 done"
} > $L 2>&1
