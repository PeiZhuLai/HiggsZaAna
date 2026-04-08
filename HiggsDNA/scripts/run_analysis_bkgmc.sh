#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# HiggsDNA will skip sample/year pairs that do not exist in the catalog,
# so we can safely submit the full Run3 DYG union in one go.
export SAMPLE_LIST="${SAMPLE_LIST:-DYGto2LG_10to50,DYGto2LG_50to100,DYJetsToLL,DYGto2LG_10to100,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau}"
export YEARS="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
export SHORT="${SHORT:-0}"

exec "${script_dir}/run_analysis_bkgmc_tmp.sh" "$@"
