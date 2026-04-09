#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}/.."

run_signal="${RUN_SIGNAL:-1}"
run_bkg="${RUN_BKG:-1}"
run_data="${RUN_DATA:-1}"

if [[ "${run_signal}" == "1" ]]; then
    echo "[run_ana_sigbkgdata] Submitting signal MC"
    bash "${script_dir}/run_analysis_signal.sh" "$@"
fi

if [[ "${run_bkg}" == "1" ]]; then
    echo "[run_ana_sigbkgdata] Submitting background MC"
    bash "${script_dir}/run_analysis_bkgmc.sh" "$@"
fi

if [[ "${run_data}" == "1" ]]; then
    echo "[run_ana_sigbkgdata] Submitting data"
    bash "${script_dir}/run_analysis_data.sh" "$@"
fi
