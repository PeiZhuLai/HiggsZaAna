#!/usr/bin/env bash
# Converts only the two DYGto2LG_10to100 2022 eras that the p2root run could not do,
# because 1_run_P2Root.sh still asked for the retired PTG-10to50 / PTG-50to100 slices.
# Everything else from that run (1205 files) already passed GATE 3.
set -uo pipefail
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
IN=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC
OUT=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix/DYGto2LG_10to100
mkdir -p "$OUT"
for era in 2022preEE 2022postEE; do
  echo "=== $era ==="
  python /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py \
    -i "$IN/DYGto2LG_10to100_${era}/merged_nominal.parquet" \
    -o "$OUT/${era}.root" --split || echo "FAILED $era"
done
echo "=== products ==="
ls -la "$OUT"
