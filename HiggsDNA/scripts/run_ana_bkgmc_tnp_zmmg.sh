#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

# Keep ROOT/cling isolated from user-site Python packages and stray include paths.
export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset CPATH
unset CPLUS_INCLUDE_PATH
unset C_INCLUDE_PATH
export PYTHONPATH="${repo_dir}"

outdir="${OUTDIR:-/eos/home-p/pelai/HZa/parquet_tnp_zmmg_tmp/mc}"
config="${CONFIG:-metadata/za_mc_tnp_zmmg_run3.json}"
log_level="${LOG_LEVEL:-INFO}"
n_cores="${N_CORES:-10}"
batch_system="${BATCH_SYSTEM:-condor}"
# batch_system="${BATCH_SYSTEM:-local}"
sample_list="${SAMPLE_LIST:-DYJetsToLL,DYJetsTo2Mu,DYJetsToLL_MLM}"
years="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
clean_analysis_state="${CLEAN_ANALYSIS_STATE:-1}"
unretire_jobs="${UNRETIRE_JOBS:-1}" # 1 re-run unfinished jobs # 0 no re-run
retire_jobs="${RETIRE_JOBS:-0}" # 1 merged parquet files # 0 no merged until all jobs finished
reconfigure_jobs="${RECONFIGURE_JOBS:-1}"
merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge
short="${SHORT:-0}"
dry_run="${DRY_RUN:-0}"

if [[ "${clean_analysis_state}" == "1" ]]; then
    rm -f "${outdir}/analysis_manager.pkl" "${outdir}/analysis_manager_temp.pkl"
fi

cmd=(
    python -s
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

if [[ "${reconfigure_jobs}" == "1" ]]; then
    cmd+=(--reconfigure_jobs)
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

cmd+=("$@")

if [[ "${dry_run}" == "1" ]]; then
    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    printf '\n'
    exit 0
fi

"${cmd[@]}"
