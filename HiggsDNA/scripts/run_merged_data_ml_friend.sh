#!/usr/bin/env bash
# Run HiggsDNA merged-ML tagger with the MLPhoton friend-tree on 2024 DATA,
# one per-dataset tag at a time. Mirrors run_merged_friend_all.sh but points at
# the per-dataset data catalog + the data-ML config (empty systematics, is_data).
#
# Each tag == one central NanoAODv15 dataset (e.g. Data_2024_EGamma0_RunC_v1),
# its MLNanoAOD friend dir, and its parent map. Produces a friend parquet with
# MLPhoton_lead_* + pass_allcuts_merged_ML so downstream can apply the per-mA
# ROI on MLPhoton_lead_mass and fit H_merged_ML_mass (= Z + MLPhoton_lead).
#
# Usage:
#   bash scripts/run_merged_data_ml_friend.sh                       # all 32 tags
#   TAGS=Data_2024_EGamma0_RunC_v1 SHORT=1 BATCH_SYSTEM=local \
#        bash scripts/run_merged_data_ml_friend.sh                  # smoke test
#
# Environment:
#   TAGS            comma-separated tag list (default = all 32 in the catalog)
#   BATCH_SYSTEM    condor (default) / local
#   N_CORES         per-sample cores (default 4)
#   FPO             files per job (default 2; data -> keep LOW, see progress doc s.10)
#   SHORT           1 to enable SHORT mode (one file per sample, fast test)
#   OUTDIR_BASE     base output dir (default /eos/.../parquet_friend)
#   BASE_CONFIG     data-ML config to patch (default metadata/za_merged_data_ML_2024.json)
#   CONDOR_REQ_MEMORY  MB for Condor RequestMemory (default 8000)

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

unset PYTHONPATH
export PYTHONPATH="${repo_dir}"

CATALOG="${CATALOG:-metadata/samples/za_merged_data_2024_perds.json}"
DEFAULT_TAGS="$(python3 -c "import json;print(','.join(sorted(json.load(open('${CATALOG}')).keys())))")"

TAGS="${TAGS:-${DEFAULT_TAGS}}"
BATCH_SYSTEM="${BATCH_SYSTEM:-condor}"
N_CORES="${N_CORES:-4}"
FPO="${FPO:-2}"
SHORT="${SHORT:-0}"
OUTDIR_BASE="${OUTDIR_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend}"
BASE_CONFIG="${BASE_CONFIG:-metadata/za_merged_data_ML_2024.json}"
CONDOR_REQ_MEMORY="${CONDOR_REQ_MEMORY:-8000}"
MLNANO_BASE="/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD"

export HIGGSDNA_CONDOR_REQ_MEMORY="${CONDOR_REQ_MEMORY}"
export HIGGSDNA_PARQUET_READ_RETRIES="${HIGGSDNA_PARQUET_READ_RETRIES:-6}"
export HIGGSDNA_PARQUET_READ_RETRY_DELAY="${HIGGSDNA_PARQUET_READ_RETRY_DELAY:-20}"
# B-mode: ship the conda-pack tarball from eos-cms and activate it node-locally on
# each worker (avoids the eos-home fuse I/O storm). getenv=True in the condor submit
# template propagates this to the job; jobs.py sees it and skips the eos-path python
# rewrite. Set to empty to fall back to the classic eos-fuse env.
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"

CFG_DIR_REL="metadata/friend_tmp_configs"
mkdir -p "${repo_dir}/${CFG_DIR_REL}" "${OUTDIR_BASE}"

n_done=0; n_fail=0
IFS=',' read -ra TAG_LIST <<< "${TAGS}"
for tag in "${TAG_LIST[@]}"; do
    [ -z "${tag}" ] && continue
    parent_map="metadata/parent_maps/${tag}.json"
    ml_dir="${MLNANO_BASE}/${tag}"
    outdir="${OUTDIR_BASE}/${tag}"

    # Resumability: skip a tag already merged (so a re-launch after an interruption
    # only redoes unfinished tags instead of the whole 36-tag sweep). Set
    # SKIP_MERGED=0 to force reprocessing.
    if [ "${SKIP_MERGED:-1}" = "1" ] && [ -f "${outdir}/merged_nominal.parquet" ]; then
        echo "[skip] ${tag} already merged (${outdir}/merged_nominal.parquet)"
        n_done=$((n_done + 1)); continue
    fi

    if [ ! -f "${parent_map}" ]; then
        echo "[skip] missing parent map: ${parent_map}  (build with scripts/build_parent_map.py)"
        n_fail=$((n_fail + 1)); continue
    fi
    if [ ! -d "${ml_dir}" ]; then
        echo "[skip] missing MLNanoAOD dir: ${ml_dir}"
        n_fail=$((n_fail + 1)); continue
    fi

    cfg_rel="${CFG_DIR_REL}/${tag}.json"
    cfg="${repo_dir}/${cfg_rel}"
    python3 - "${BASE_CONFIG}" "${cfg}" "${tag}" "${parent_map}" "${ml_dir}" "${CATALOG}" <<'PY'
import json, os, sys
base, out, tag, parent_map, ml_dir, catalog = sys.argv[1:]
with open(base) as f:
    d = json.load(f)
# era-generic: read the tag's era key from the per-ds catalog (2024 -> "2024",
# 2223 -> "2022preEE"/... ); also pin samples.catalog so the base config's
# catalog field cannot mismatch the tag's era (works for 2024 and 2022/2023).
cat_d = json.load(open(catalog))
era = list(cat_d[tag]["files"].keys())[0]
d["samples"]["catalog"] = catalog
d["samples"]["sample_list"] = [tag]
d["samples"]["years"] = [era]
d["mlphoton_friend"] = True
# ABSOLUTE path so condor workers (cwd = job scratch dir) can find the parent map on afs;
# the lxplus submit template only ships GRID_PROXY, executable+config are read from afs.
d["mlphoton_parent_map"] = os.path.abspath(parent_map)
d["mlphoton_dir"] = ml_dir
d["mlphoton_tag"] = tag
with open(out, "w") as f:
    json.dump(d, f, indent=2)
PY

    echo "=== ${tag} (fpo=${FPO}) ==="
    mkdir -p "${outdir}"
    rm -f "${outdir}/analysis_manager.pkl" "${outdir}/analysis_manager_temp.pkl"

    cmd=(
        python scripts/run_analysis.py
        --config "${cfg_rel}"
        --log-level INFO
        --n_cores "${N_CORES}"
        --output_dir "${outdir}"
        --batch_system "${BATCH_SYSTEM}"
        --merge_outputs
        --fpo "${FPO}"
    )
    [ "${SHORT}" = "1" ] && cmd+=(--short)

    if "${cmd[@]}"; then
        n_done=$((n_done + 1))
    else
        echo "[FAIL] ${tag}"; n_fail=$((n_fail + 1))
    fi
done

echo
echo "Done: ${n_done} tags; failed: ${n_fail}"
