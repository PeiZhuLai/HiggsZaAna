#!/usr/bin/env bash
# Produce the 2024 merged-selection DATA friend parquet -- the missing input
# that unblocks the data-driven background fit in ../1_grand_merged.sh.
#
# DATA does NOT need the MLPhoton regression (the merged selection
# pass_allcuts_merged_AN2020 and H_mass come from central NanoAODv15), so this
# runs the merged tagger in plain (non-friend) mode on the 2024 EGamma datasets
# already catalogued under "Data"/2024 in metadata/samples/zgamma_tutorial.json.
#
# It delegates to run_merged_data.sh (which carries the resubmit/merge logic)
# but pins the 2024-only, ML-stripped config and an output path that
# prep_bkg_data.py auto-detects (parquet_friend/Data_2024/...).
#
# Usage:
#   bash scripts/run_merged_data_2024.sh                 # condor, full 2024
#   SHORT=1 bash scripts/run_merged_data_2024.sh         # 1 file, smoke test
#   BATCH_SYSTEM=local bash scripts/run_merged_data_2024.sh
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

# higgs_dna is imported from the repo (not pip-installed in hza_ana); the plain
# data runner does not set this, so the wrapper must (mirrors the friend driver).
export PYTHONPATH="${repo_dir}:${PYTHONPATH:-}"

export CONFIG="${CONFIG:-metadata/za_merged_data_2024.json}"
export OUTDIR="${OUTDIR:-/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend/Data_2024}"
export YEARS="${YEARS:-2024}"
export SAMPLE_LIST="${SAMPLE_LIST:-Data}"
export BATCH_SYSTEM="${BATCH_SYSTEM:-condor}"
export N_CORES="${N_CORES:-4}"
export FPO="${FPO:-2}"
export MERGE_OUTPUTS="${MERGE_OUTPUTS:-1}"
export UNRETIRE_JOBS="${UNRETIRE_JOBS:-1}"
export SHORT="${SHORT:-0}"
# condor hardening (ref: doc §8.2b / HZa condor hardening): fresh jobs only
export CONDOR_REQ_MEMORY="${CONDOR_REQ_MEMORY:-13000}"

echo "[run_merged_data_2024] config=${CONFIG}"
echo "[run_merged_data_2024] outdir=${OUTDIR}"
echo "[run_merged_data_2024] years=${YEARS}  batch=${BATCH_SYSTEM}  short=${SHORT}"

exec bash scripts/run_merged_data.sh "$@"
