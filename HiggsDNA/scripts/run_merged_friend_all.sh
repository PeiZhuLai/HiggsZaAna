#!/usr/bin/env bash
# Run HiggsDNA merged tagger with MLPhoton friend-tree for all 19 signal +
# 4 bkg tags, one sample at a time. Each invocation rewrites the top-level
# mlphoton_* keys in a temp config and calls run_analysis.py.
#
# Usage:
#   bash scripts/run_merged_friend_all.sh                # all 23 tags
#   TAGS=mA_M0p4,mA_M5 bash scripts/run_merged_friend_all.sh
#   BATCH_SYSTEM=condor bash scripts/run_merged_friend_all.sh
#   OUTDIR=/eos/.../parquet_friend bash scripts/run_merged_friend_all.sh
#
# Environment:
#   TAGS            comma-separated tag list (default = signal + 2024 bkg)
#   BATCH_SYSTEM    condor (default) / local
#   N_CORES         per-sample cores (default 4)
#   SHORT           1 to enable SHORT mode (one file per sample, fast test)
#   OUTDIR_BASE     base output dir (default /eos/.../parquet_friend)
#   BASE_CONFIG     base config to patch (default metadata/za_merged_signal_run3.json)

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

unset PYTHONPATH
export PYTHONPATH="${repo_dir}"

DEFAULT_TAGS="mA_M0p1,mA_M0p2,mA_M0p3,mA_M0p4,mA_M0p5,mA_M0p6,mA_M0p7,mA_M0p8,mA_M0p9,mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,Bkg_DYGto2LG_10to100_2024,Bkg_DYJetsTo2E_2024,Bkg_DYJetsTo2Mu_2024,Bkg_DYJetsTo2Tau_2024"

TAGS="${TAGS:-${DEFAULT_TAGS}}"
BATCH_SYSTEM="${BATCH_SYSTEM:-condor}"
N_CORES="${N_CORES:-4}"
SHORT="${SHORT:-0}"
OUTDIR_BASE="${OUTDIR_BASE:-/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_friend}"
BASE_CONFIG="${BASE_CONFIG:-metadata/za_merged_signal_run3.json}"

# Map each tag to (sample_tag in zgamma_tutorial.json, mlphoton_dir)
# Signal: sample_tag = tag, mlphoton_dir = /eos/.../MLNanoAOD/<tag>
# Bkg: sample_tag = tag without Bkg_ prefix and _2024 suffix (e.g. DYJetsTo2E)
declare -A SAMPLE_TAG_MAP
SAMPLE_TAG_MAP[Bkg_DYGto2LG_10to100_2024]="DYGto2LG_10to100"
SAMPLE_TAG_MAP[Bkg_DYJetsTo2E_2024]="DYJetsTo2E"
SAMPLE_TAG_MAP[Bkg_DYJetsTo2Mu_2024]="DYJetsTo2Mu"
SAMPLE_TAG_MAP[Bkg_DYJetsTo2Tau_2024]="DYJetsTo2Tau"

mkdir -p "${OUTDIR_BASE}"
# Configs must live under the HiggsDNA repo because run_analysis.py's
# expand_path() blindly prepends <HiggsDNA root>/ to whatever is passed via
# --config, so an absolute /tmp path is mangled into <repo>//tmp/...
CFG_DIR_REL="metadata/friend_tmp_configs"
mkdir -p "${repo_dir}/${CFG_DIR_REL}"

n_done=0; n_fail=0
IFS=',' read -ra TAG_LIST <<< "${TAGS}"
for tag in "${TAG_LIST[@]}"; do
    [ -z "${tag}" ] && continue
    sample_tag="${SAMPLE_TAG_MAP[${tag}]:-${tag}}"
    parent_map="metadata/parent_maps/${tag}.json"
    ml_dir="/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/${tag}"
    outdir="${OUTDIR_BASE}/${tag}"

    if [ ! -f "${parent_map}" ]; then
        echo "[skip] missing parent map: ${parent_map}"
        n_fail=$((n_fail + 1))
        continue
    fi
    if [ ! -d "${ml_dir}" ]; then
        echo "[skip] missing MLNanoAOD dir: ${ml_dir}"
        n_fail=$((n_fail + 1))
        continue
    fi

    # Patch base config: set sample_list, mlphoton_*, then write to a temp file
    # (path is repo-relative so run_analysis.py's expand_path resolves it correctly)
    cfg_rel="${CFG_DIR_REL}/${tag}.json"
    cfg="${repo_dir}/${cfg_rel}"
    python3 - "${BASE_CONFIG}" "${cfg}" "${sample_tag}" "${parent_map}" "${ml_dir}" "${tag}" <<'PY'
import json, sys
base, out, sample_tag, parent_map, ml_dir, friend_tag = sys.argv[1:]
with open(base) as f:
    d = json.load(f)
d["samples"]["sample_list"] = [sample_tag]
d["samples"]["years"] = ["2024"]
d["mlphoton_friend"] = True
d["mlphoton_parent_map"] = parent_map
d["mlphoton_dir"] = ml_dir
d["mlphoton_tag"] = friend_tag
with open(out, "w") as f:
    json.dump(d, f, indent=2)
PY

    echo "=== ${tag} (sample_tag=${sample_tag}) ==="
    mkdir -p "${outdir}"
    rm -f "${outdir}/analysis_manager.pkl" "${outdir}/analysis_manager_temp.pkl"

    cmd=(
        python scripts/run_analysis.py
        --config "${cfg_rel}"
        --log-level INFO
        --n_cores "${N_CORES}"
        --output_dir "${outdir}"
        --batch_system "${BATCH_SYSTEM}"
        --unretire_jobs
        --merge_outputs
        --fpo 1
    )
    [ "${SHORT}" = "1" ] && cmd+=(--short)

    if "${cmd[@]}"; then
        n_done=$((n_done + 1))
    else
        echo "[FAIL] ${tag}"
        n_fail=$((n_fail + 1))
    fi
done

echo
echo "Done: ${n_done} tags; failed: ${n_fail}"
