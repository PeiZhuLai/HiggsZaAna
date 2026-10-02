#!/bin/bash
# S1 (L3 review, Si Hyun Jeon): derive alternative sideband reweights.
# Same inputs/options as the nominal FSR-fix derivation (all script defaults:
# run3_bdt_inputs_fsrfix Data/All_Bkg run3.root, tree inclusive, 5 iterations,
# 30 quantile bins, clip [0.2,5], seed 12345, weight 'factor'); only --sideband,
# the output JSON and the plot dir differ. Never writes the nominal JSON.
# Usage: run_derive.sh <sideband> <out_json> [extra args]
set -euo pipefail
SB=$1; OUT=$2; shift 2
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
SCR=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/1_make_sideband_reweight.py
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927
case "$OUT" in *sideband_run3_iterative_fsrfix.json|*sideband_run3_iterative.json) echo "refuse to overwrite nominal"; exit 3;; esac
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
env -i HOME=$HOME PATH=/usr/bin:/bin $PY $SCR --sideband "$SB" -o "$OUT" --plot-dir "$S/validation_pdfs/$SB" "$@"
echo "DONE $SB $(date)"
