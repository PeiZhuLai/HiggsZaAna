#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$script_dir/make_dataVmc_condor_jobs.py" --condor-dir "$script_dir" "$@"
condor_submit "$script_dir/dataVmc.submit"
