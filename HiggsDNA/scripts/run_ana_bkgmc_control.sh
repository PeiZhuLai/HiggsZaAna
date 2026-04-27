#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

outdir="${OUTDIR:-/eos/home-p/pelai/HZa/parquet_DNA_control_tmp/Bkg_MC_control}"
config="${CONFIG:-metadata/za_bkgmc_run3_control.json}"
log_level="${LOG_LEVEL:-INFO}" # DEBUG # INFO
n_cores="${N_CORES:-10}"
batch_system="${BATCH_SYSTEM:-condor}" # local # condor
sample_list="${SAMPLE_LIST:-DYGto2LG_10to50,DYGto2LG_50to100,DYJetsToLL,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau}"
years="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
# sample_list="${SAMPLE_LIST:-DYGto2LG_10to50,DYGto2LG_50to100,DYJetsToLL,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau}"
# years="${YEARS:-2024}"
fpo="${FPO:-5}"
clean_analysis_state="${CLEAN_ANALYSIS_STATE:-0}"
unretire_jobs="${UNRETIRE_JOBS:-1}"
merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge
short="${SHORT:-0}" # 1 short # 0 full
dry_run="${DRY_RUN:-0}"
condor_submit_chunk_size="${CONDOR_SUBMIT_CHUNK_SIZE:-500}"
condor_q_timeout="${CONDOR_Q_TIMEOUT:-60}"

export HIGGSDNA_CONDOR_SUBMIT_CHUNK_SIZE="${condor_submit_chunk_size}"
export HIGGSDNA_CONDOR_Q_TIMEOUT="${condor_q_timeout}"

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

if [[ "${merge_outputs}" == "1" ]]; then
    cmd+=(--merge_outputs)
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

if [[ -n "${fpo}" ]]; then
    cmd+=(--fpo "${fpo}")
fi

cmd+=("$@")

if [[ "${dry_run}" == "1" ]]; then
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

"${cmd[@]}"
