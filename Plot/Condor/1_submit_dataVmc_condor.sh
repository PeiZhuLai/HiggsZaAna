#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${SKIP_ENV_PACK:-0}" != "1" ]]; then
    "$script_dir/pack_conda_env_for_condor.sh"
fi

python3 "$script_dir/make_dataVmc_condor_jobs.py" --condor-dir "$script_dir" "$@"
condor_submit "$script_dir/dataVmc.submit"
