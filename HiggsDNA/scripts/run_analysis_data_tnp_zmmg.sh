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

outdir="/eos/home-p/pelai/HZa/parquet_tnp_zmmg/data"
reconfigure_jobs="${RECONFIGURE_JOBS:-1}"

# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/data/job_1
# rm -fr /eos/home-p/pelai/HZa/parquet_tnp_zmmg/data/analysis_manager.pkl

cmd=(
    python -s
    scripts/run_analysis.py
    --config "metadata/za_data_tnp_zmmg_run3.json"
    --log-level "INFO"
    --n_cores 10
    --output_dir "$outdir"
    --unretire_jobs
    --batch_system "condor"
)

if [[ "${reconfigure_jobs}" == "1" ]]; then
    cmd+=(--reconfigure_jobs)
fi

"${cmd[@]}" #--short #--batch_system "local" "condor" INFO DEBUG #--with_skimmed
