#!/usr/bin/env bash
# 2026-10-01: rescore the V3/V8 extra-background samples with the BDT retrained after the 2024
# DY+jets overlap veto, and recompute the yield table against the new (vetoed) DY+jets.
# Old scored files kept in run3_bdt_scored_fsrfix_extraBkg_preDYveto_20261001; old tables as *.preDYveto_20261001.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927
P=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkg2024
R=/eos/home-p/pelai/HZa/root_P2Root
SC=$R/run3_bdt_scored_fsrfix_extraBkg
[ -d $R/run3_bdt_scored_fsrfix_extraBkg_preDYveto_20261001 ] || mv $SC $R/run3_bdt_scored_fsrfix_extraBkg_preDYveto_20261001
for f in yields_extrabkg.md yields_extrabkg.json; do [ -f $D/$f ] && cp -n $D/$f $D/$f.preDYveto_20261001; done
for s in TTto2L2Nu TTG TTGG_Run3 ZGG; do
  bash $D/score_local.sh $P/${s}_2024/merged_nominal.parquet $SC/$s/2024.root --sideband-reweight-mode always > $D/logs/rescore_dyveto_$s.log 2>&1
  echo "$s rc=$? $(ls -la $SC/$s/2024.root 2>&1 | awk '{print $5}')"
done
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
env -i HOME=$HOME PATH=/usr/bin:/bin $PY $D/compute_extrabkg_yields.py --new-base $SC --out-prefix $D/yields_extrabkg > $D/logs/yields_dyveto.log 2>&1
echo "yields rc=$?"
