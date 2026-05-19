#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

outdir="${OUTDIR:-/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Data}"
config="${CONFIG:-metadata/za_data_run3.json}"
log_level="${LOG_LEVEL:-INFO}" # DEBUG # INFO
n_cores="${N_CORES:-10}"
batch_system="${BATCH_SYSTEM:-condor}" # local # condor
sample_list="${SAMPLE_LIST:-Data}"
years="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
fpo="${FPO:-4}"
clean_analysis_state="${CLEAN_ANALYSIS_STATE:-1}"
unretire_jobs="${UNRETIRE_JOBS:-1}" # 1 re-run unfinished jobs # 0 no re-run
retire_jobs="${RETIRE_JOBS:-0}" # 1 merged parquet files # 0 no merged until all jobs finished
merge_outputs="${MERGE_OUTPUTS:-1}" # 1 merge # 0 no merge
short="${SHORT:-0}" # 1 short # 0 full
dry_run="${DRY_RUN:-0}"
reconfigure_jobs="${RECONFIGURE_JOBS:-0}" # 1 rewrite job configs/scripts/condor submit files
condor_req_memory="${CONDOR_REQ_MEMORY:-20000}" # MB, passed to Condor RequestMemory
parquet_read_retries="${HIGGSDNA_PARQUET_READ_RETRIES:-6}"
parquet_read_retry_delay="${HIGGSDNA_PARQUET_READ_RETRY_DELAY:-10}"

export HIGGSDNA_CONDOR_REQ_MEMORY="${condor_req_memory}"
export HIGGSDNA_PARQUET_READ_RETRIES="${parquet_read_retries}"
export HIGGSDNA_PARQUET_READ_RETRY_DELAY="${parquet_read_retry_delay}"

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

if [[ "${retire_jobs}" == "1" ]]; then
    cmd+=(--retire_jobs)
elif [[ "${unretire_jobs}" == "1" ]]; then
    cmd+=(--unretire_jobs)
fi

if [[ "${merge_outputs}" == "1" ]]; then
    cmd+=(--merge_outputs)
fi

if [[ "${short}" == "1" ]]; then
    cmd+=(--short)
fi

if [[ "${reconfigure_jobs}" == "1" ]]; then
    cmd+=(--reconfigure_jobs)
fi

if [[ -n "${sample_list}" ]]; then
    cmd+=(--sample_list "${sample_list}")
fi

if [[ -n "${years}" ]]; then
    cmd+=(--years "${years}")
fi

if [[ -n "${fpo}" ]]; then
    cmd+=(--fpo "${fpo}")
fi

cmd+=("$@")

if [[ "${dry_run}" == "1" ]]; then
    printf 'HIGGSDNA_CONDOR_REQ_MEMORY=%q\n' "${HIGGSDNA_CONDOR_REQ_MEMORY}"
    printf 'HIGGSDNA_PARQUET_READ_RETRIES=%q\n' "${HIGGSDNA_PARQUET_READ_RETRIES}"
    printf 'HIGGSDNA_PARQUET_READ_RETRY_DELAY=%q\n' "${HIGGSDNA_PARQUET_READ_RETRY_DELAY}"
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

"${cmd[@]}"
