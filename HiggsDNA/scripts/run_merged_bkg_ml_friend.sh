#!/usr/bin/env bash
# Run the HiggsDNA merged-ML tagger WITH the MLPhoton friend tree on 2024 background MC.
#
# Why this exists: run_merged_bkgmc.sh drives the same tagger but its config
# (metadata/za_merged_bkgmc_run3.json) never sets mlphoton_friend, so the join is
# simply never switched on. The resulting friend parquet has 263 columns and no
# MLPhoton_* at all -- only the older pass_allcuts_merged_AN2020. Nothing errors:
# add_merged_flag.py fills the missing pass_allcuts_merged_ML with 0 and the log
# reads "pass_merged 0", which looks exactly like "no background event passes the
# merged selection". It is not -- the information was never produced.
#
# This mirrors run_merged_data_ml_friend.sh: one tag per (sample, era), each with
# its own MLNanoAOD dir and parent map, patched into a per-tag config.
#
# Tag convention: Bkg_<sample>_<era>, matching MLNanoAOD/ and metadata/parent_maps/.
#
# Usage:
#   bash scripts/run_merged_bkg_ml_friend.sh                          # all 4 tags
#   TAGS=Bkg_DYJetsTo2E_2024 SHORT=1 BATCH_SYSTEM=local \
#        bash scripts/run_merged_bkg_ml_friend.sh                     # smoke test
#   DRY_RUN=1 bash scripts/run_merged_bkg_ml_friend.sh                # print, submit nothing
#
# Environment:
#   TAGS            comma-separated (default: the four 2024 backgrounds)
#   BATCH_SYSTEM    condor (default) / local
#   N_CORES         per-sample cores (default 4)
#   FPO             files per job (default 4)
#   SHORT           1 = one file per sample
#   DRY_RUN         1 = show what would run, submit nothing
#   OUTDIR_BASE     default /eos/.../HZa_merged/parquet_friend
#   BASE_CONFIG     MC config to patch (default metadata/za_merged_bkgmc_run3.json)
#   SKIP_MERGED     1 (default) = skip a tag whose merged_nominal.parquet exists

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

unset PYTHONPATH
export PYTHONPATH="${repo_dir}"

TAGS="${TAGS:-Bkg_DYGto2LG_10to100_2024,Bkg_DYJetsTo2E_2024,Bkg_DYJetsTo2Mu_2024,Bkg_DYJetsTo2Tau_2024}"
BATCH_SYSTEM="${BATCH_SYSTEM:-condor}"
N_CORES="${N_CORES:-4}"
FPO="${FPO:-4}"
SHORT="${SHORT:-0}"
DRY_RUN="${DRY_RUN:-0}"
OUTDIR_BASE="${OUTDIR_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML}"
BASE_CONFIG="${BASE_CONFIG:-metadata/za_merged_bkgmc_run3.json}"
CONDOR_REQ_MEMORY="${CONDOR_REQ_MEMORY:-8000}"
MLNANO_BASE="${MLNANO_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD}"

export HIGGSDNA_CONDOR_REQ_MEMORY="${CONDOR_REQ_MEMORY}"
export HIGGSDNA_PARQUET_READ_RETRIES="${HIGGSDNA_PARQUET_READ_RETRIES:-6}"
export HIGGSDNA_PARQUET_READ_RETRY_DELAY="${HIGGSDNA_PARQUET_READ_RETRY_DELAY:-20}"
# B-mode conda-pack, same as the data friend run -- avoids the eos-home fuse storm.
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"

CFG_DIR_REL="metadata/friend_tmp_configs"
mkdir -p "${repo_dir}/${CFG_DIR_REL}" "${OUTDIR_BASE}"

n_done=0; n_fail=0; n_skip=0
IFS=',' read -ra TAG_LIST <<< "${TAGS}"
for tag in "${TAG_LIST[@]}"; do
    [ -z "${tag}" ] && continue

    # Bkg_<sample>_<era> -> sample, era
    era="${tag##*_}"
    sample="${tag#Bkg_}"; sample="${sample%_${era}}"

    parent_map="metadata/parent_maps/${tag}.json"
    ml_dir="${MLNANO_BASE}/${tag}"
    outdir="${OUTDIR_BASE}/${tag}"

    if [ "${SKIP_MERGED:-1}" = "1" ] && [ -f "${outdir}/merged_nominal.parquet" ]; then
        echo "[skip] ${tag} already merged (${outdir}/merged_nominal.parquet)"
        n_skip=$((n_skip + 1)); continue
    fi
    if [ ! -f "${parent_map}" ]; then
        echo "[skip] missing parent map: ${parent_map}"
        n_fail=$((n_fail + 1)); continue
    fi
    if [ ! -d "${ml_dir}" ]; then
        echo "[skip] missing MLNanoAOD dir: ${ml_dir}"
        n_fail=$((n_fail + 1)); continue
    fi

    # An all-empty parent map is the failure mode this whole script exists to avoid
    # (DAS returns no file-level parentage for some datasets -- see
    # build_parent_map_lumi.py). Refuse rather than produce MLPhoton-less parquet.
    n_links=$(python3 -c "import json;d=json.load(open('${parent_map}'));print(sum(len(v) for v in d.values()))")
    if [ "${n_links}" = "0" ]; then
        echo "[FAIL] ${tag}: parent map has 0 links -- rebuild with scripts/build_parent_map_lumi.py"
        n_fail=$((n_fail + 1)); continue
    fi

    cfg_rel="${CFG_DIR_REL}/${tag}.json"
    cfg="${repo_dir}/${cfg_rel}"
    python3 - "${BASE_CONFIG}" "${cfg}" "${tag}" "${parent_map}" "${ml_dir}" "${sample}" "${era}" <<'PY'
import json, os, sys
base, out, tag, parent_map, ml_dir, sample, era = sys.argv[1:]
with open(base) as f:
    d = json.load(f)
d["samples"]["sample_list"] = [sample]
d["samples"]["years"] = [era]
d["mlphoton_friend"] = True
# ABSOLUTE -- a relative path leaves the loader unable to find the map on a condor
# worker, and it then skips the join silently instead of failing.
d["mlphoton_parent_map"] = os.path.abspath(parent_map)
d["mlphoton_dir"] = ml_dir
d["mlphoton_tag"] = tag
with open(out, "w") as f:
    json.dump(d, f, indent=2)
print(f"  config: {out}  (sample={sample}, era={era}, links ok)")
PY

    echo "=== ${tag}  sample=${sample} era=${era}  fpo=${FPO} links=${n_links} ==="
    if [ "${DRY_RUN}" = "1" ]; then
        echo "  would run: python scripts/run_analysis.py --config ${cfg_rel} --output_dir ${outdir} --batch_system ${BATCH_SYSTEM} --fpo ${FPO}"
        n_done=$((n_done + 1)); continue
    fi

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
echo "Done: ${n_done}; skipped: ${n_skip}; failed: ${n_fail}"
