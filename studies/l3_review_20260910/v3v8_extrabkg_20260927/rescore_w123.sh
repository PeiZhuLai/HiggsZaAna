#!/usr/bin/env bash
# Rescore the V3/V8 extra-background samples (tt, ttG, ttGG, ZGG; 2024) with the current production models
# (low-mass w123, high-mass unchanged) and recompute the yield table at the current working points.
# Input : parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkg2024/<s>_2024/merged_nominal.parquet
# Output: run3_bdt_scored_fsrfix_extraBkg_w123/<s>/2024.root (the 10-01 files in run3_bdt_scored_fsrfix_extraBkg are kept);
#         yields_extrabkg_w123.{json,md}; logs/rescore_w123_<s>.log, logs/yields_w123.log
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927
P=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkg2024
SC=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_extraBkg_w123
for s in TTto2L2Nu TTG TTGG_Run3 ZGG; do
  bash $D/score_local.sh $P/${s}_2024/merged_nominal.parquet $SC/$s/2024.root --sideband-reweight-mode always > $D/logs/rescore_w123_$s.log 2>&1
  rc=$?; echo "$s rc=$rc $(ls -la $SC/$s/2024.root 2>&1 | awk '{print $5}')"
  [ $rc -eq 0 ] || { echo FAILED; exit 1; }
done
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
env -i HOME=$HOME PATH=/usr/bin:/bin $PY $D/compute_extrabkg_yields.py --new-base $SC --out-prefix $D/yields_extrabkg_w123 > $D/logs/yields_w123.log 2>&1
echo "yields rc=$?"
