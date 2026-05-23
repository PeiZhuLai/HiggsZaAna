#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${SKIP_ENV_PACK:-0}" != "1" ]]; then
    "$script_dir/pack_conda_env_for_condor.sh"
fi

has_final_tags=0
for arg in "$@"; do
    case "$arg" in
        --final-tags|--final-tags=*) has_final_tags=1 ;;
    esac
done

make_args=(--condor-dir "$script_dir")
if [[ "$has_final_tags" == "0" ]]; then
    make_args+=(--final-tags sideband_rwgt)
fi

python3 "$script_dir/make_dataVmc_condor_jobs.py" "${make_args[@]}" "$@"
if [[ "${NO_SUBMIT:-0}" == "1" ]]; then
    echo "[Info] NO_SUBMIT=1; generated dataVmc.submit and dataVmc_jobs.txt only."
    exit 0
fi
condor_submit "$script_dir/dataVmc.submit"
