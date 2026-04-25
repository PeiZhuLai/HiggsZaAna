#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

outdir="${OUTDIR:-/eos/home-p/pelai/HZa/parquet_DNA_tmp/Sig_MC}"
config="${CONFIG:-metadata/za_signal_run3.json}"
log_level="${LOG_LEVEL:-DEBUG}"
n_cores="${N_CORES:-15}"
batch_system="${BATCH_SYSTEM:-local}"
# sample_list="${SAMPLE_LIST:-mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30}"
# years="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
sample_list="${SAMPLE_LIST:-mA_M1}"
years="${YEARS:-2024}"
clean_analysis_state="${CLEAN_ANALYSIS_STATE:-1}"
unretire_jobs="${UNRETIRE_JOBS:-1}"
short="${SHORT:-1}"
dry_run="${DRY_RUN:-0}"

if [[ "${clean_analysis_state}" == "1" ]]; then
    rm -f "${outdir}/analysis_manager.pkl" "${outdir}/analysis_manager_temp.pkl"
fi

cmd=(
    python
    scripts/run_analysis.py
    --config "${config}"
    --log-level "${log_level}"
    --n_cores "${n_cores}"
    --output_dir "${outdir}"
    --batch_system "${batch_system}"
)

if [[ "${unretire_jobs}" == "1" ]]; then
    cmd+=(--unretire_jobs)
fi

if [[ "${short}" == "1" ]]; then
    cmd+=(--short)
fi

if [[ -n "${sample_list}" ]]; then
    cmd+=(--sample_list "${sample_list}")
fi

if [[ -n "${years}" ]]; then
    cmd+=(--years "${years}")
fi

cmd+=("$@")

if [[ "${dry_run}" == "1" ]]; then
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

"${cmd[@]}"
