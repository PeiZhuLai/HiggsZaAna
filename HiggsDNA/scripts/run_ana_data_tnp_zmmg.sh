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

outdir="/eos/home-p/pelai/HZa/parquet_tnp_zmmg_tmp/data"
fpo="${FPO:-5}"
unretire_jobs="${UNRETIRE_JOBS:-1}" # 1 re-run unfinished jobs # 0 no re-run
retire_jobs="${RETIRE_JOBS:-0}" # 1 merged parquet files # 0 no merged until all jobs finished
reconfigure_jobs="${RECONFIGURE_JOBS:-1}"
merge_outputs="${MERGE_OUTPUTS:-0}" # 1 merge # 0 no merge

# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg_tmp/data/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/data/analysis_manager.pkl

cmd=(
    python -s
    scripts/run_analysis.py
    --config "metadata/za_data_tnp_zmmg_run3.json"
    --log-level "INFO"
    --n_cores 10
    --output_dir "$outdir"
    --batch_system "condor"
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

if [[ -n "${fpo}" ]]; then
    cmd+=(--fpo "${fpo}")
fi

"${cmd[@]}" #--short #--batch_system "local" "condor" INFO DEBUG #--with_skimmed
