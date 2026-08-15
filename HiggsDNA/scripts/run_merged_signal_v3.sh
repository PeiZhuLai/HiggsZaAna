#!/usr/bin/env bash
# Re-run the merged signal HiggsDNA over the RETRAINED MLPhoton MLNanoAOD (v3).
#
# za_merged_signal_run3.json wires the MLPhoton friend per sample -- one
# mlphoton_dir / mlphoton_tag / mlphoton_parent_map at a time -- so a multi-mass
# production means running it once per mass point with a rewritten config. This
# script does that rewriting into a temporary config so the checked-in one is
# never edited in place (its committed values point at the OLD production and
# are what reproduce the current parquet).
#
# Usage:
#   bash run_merged_signal_v3.sh                 # the 9 sub-GeV points
#   bash run_merged_signal_v3.sh M0p1 M0p2
#   DRY_RUN=1 bash run_merged_signal_v3.sh       # print configs, submit nothing
#
# Input  : /eos/.../HZa_merged/MLNanoAOD_v3/mA_<tag>/
# Output : /eos/.../HZa_merged/parquet_merged_DNA_v3/Sig_MC_MLNANO_all/
# Configs: <repo>/metadata/_generated/za_merged_signal_v3_<tag>.json
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

# run_merged_signal.sh assumes the caller already activated an env -- it only
# fixes PYTHONPATH. `hza_ana` is the HZa HiggsDNA env (coffea/awkward/pyarrow).
# Activate non-interactively: sourcing ~/setup_*.sh relies on shell aliases and
# silently does nothing here.
if [ -z "${SKIP_ENV:-}" ]; then
    set +u
    source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
    conda activate /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
    set -u
    echo "[env] python=$(which python)"
fi

# Merged analysis only uses the sub-GeV points (1_grand_merged.sh MASSES, and
# merged_p2root.py ROI_WINDOWS both stop at 0p9).
MASSES="${*:-M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9}"

ML_BASE="${ML_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD_v3}"
OUTDIR="${OUTDIR:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_merged_DNA_v3/Sig_MC_MLNANO_all}"
BASE_CONFIG="${BASE_CONFIG:-metadata/za_merged_signal_run3.json}"
GEN_DIR="${GEN_DIR:-metadata/_generated}"
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python

mkdir -p "${GEN_DIR}" "${OUTDIR}"

for m in ${MASSES}; do
    tag="mA_${m}"
    ml_dir="${ML_BASE}/${tag}"
    pmap="metadata/parent_maps/${tag}.json"
    cfg="${GEN_DIR}/za_merged_signal_v3_${tag}.json"

    if [ ! -d "${ml_dir}" ]; then
        echo "[skip] ${tag}: no MLNanoAOD at ${ml_dir}"
        continue
    fi
    n_ml=$(find "${ml_dir}" -name '*.root' 2>/dev/null | wc -l)
    if [ ! -f "${pmap}" ]; then
        echo "[skip] ${tag}: no parent map ${pmap}"
        continue
    fi

    cat_out="${GEN_DIR}/catalog_v3_${tag}.json"
    "${PY}" - "$BASE_CONFIG" "$cfg" "$ml_dir" "$tag" "$pmap" "$cat_out" <<'PYEOF'
import json, os, sys
base, out, ml_dir, tag, pmap, cat_out = sys.argv[1:7]
d = json.load(open(base))
d["mlphoton_dir"] = os.path.abspath(ml_dir)
d["mlphoton_tag"] = tag
# ABSOLUTE. The job runs on a worker whose cwd is not the repo, so a relative
# parent-map path is simply not found -- and the loader then skips the friend
# join SILENTLY: the parquet is written, just without any MLPhoton_* column.
# (Measured: 263 columns instead of 289, no error anywhere.) run_merged_data_
# ml_friend.sh already does this; the signal path had to match.
d["mlphoton_parent_map"] = os.path.abspath(pmap)
# catalog sample names are mA_MLNANO_<tag>, the friend dirs are mA_<tag>
sample = f"mA_MLNANO_{tag.split('_', 1)[1]}"
d["samples"]["sample_list"] = [sample]
d["samples"]["years"] = ["2024"]

# The catalog must point at the NanoAODv15 DATASET, not at the MLNanoAOD.
#
# Two traps here, both of which cost a full round of failed jobs:
#   * the committed catalog points at /eos/project/.../MLNanoAOD, which stopped
#     existing when HZa_merged moved to phys_susy. An empty file list makes
#     HiggsDNA build zero jobs and die with ZeroDivisionError in the progress
#     bar -- nothing in the error mentions the missing input.
#   * MLNanoAOD is a SLIM FRIEND (18 branches: run/lumi/event + MLPhoton_*),
#     not a full NanoAOD. Pointing `files` at it makes every job fail with
#     "no field 'Photon' in record with 5 fields". The MLPhoton content reaches
#     the analysis through the friend join (mlphoton_dir + parent map), which is
#     what mlphoton_friend=True is for.
#
# The parent map's keys are the NanoAODv15 LFNs for exactly this mass point, so
# the dataset name is derived from them rather than hard-coded -- that keeps the
# nano input and the friend map guaranteed consistent.
pm_keys = list(json.load(open(pmap)).keys())
if not pm_keys:
    raise SystemExit(f"parent map {pmap} is empty")
# /store/mc/<campaign>/<primary>/<tier>/<conditions>/<block>/<file>.root
p = pm_keys[0].split("/")
campaign, primary, tier, conditions = p[3], p[4], p[5], p[6]
dataset = f"/{primary}/{campaign}-{conditions}/{tier}"

cat_path = d["samples"]["catalog"]
cat = json.load(open(cat_path))
node = cat["samples"] if "samples" in cat else cat
if sample not in node:
    raise SystemExit(f"catalog {cat_path} has no sample {sample}")
node[sample]["files"] = {era: ([dataset] if era == "2024" else [])
                         for era in node[sample]["files"]}
json.dump(cat, open(cat_out, "w"), indent=4)
# RELATIVE, unlike the parent map: SampleManager runs the catalog through
# expand_path(), which prefixes the repo dir -- an absolute path there becomes
# "<repo>//afs/...". The parent map is opened directly, so it must be absolute.
# The two fields genuinely need opposite conventions.
d["samples"]["catalog"] = cat_out
d["mlphoton_friend"] = True
print(f"    nano dataset = {dataset}")

json.dump(d, open(out, "w"), indent=4)
print(f"    sample={sample}  files -> {ml_dir}")
PYEOF

    echo "=============== ${tag}  (${n_ml} MLNanoAOD files)  -> ${cfg}"
    if [ "${DRY_RUN:-0}" = "1" ]; then
        echo "    [dry-run] would run run_merged_signal.sh with CONFIG=${cfg}"
        continue
    fi

    # One mass point must not take down the rest. HiggsDNA resolves the dataset
    # through DAS, and a transient DAS error yields an empty file list -> zero
    # jobs -> ZeroDivisionError in the progress bar. That is what killed M0p2 on
    # the first pass even though the dataset was fine (24 files on DAS), and
    # set -e then skipped every remaining mass point. Retry, then carry on.
    ok=0
    for attempt in 1 2 3; do
        set +e
        CONFIG="${cfg}" OUTDIR="${OUTDIR}" SAMPLE_LIST="mA_MLNANO_${m}" YEARS="2024" \
            bash "${script_dir}/run_merged_signal.sh"
        rc=$?
        set -e
        if [ ${rc} -eq 0 ]; then ok=1; break; fi
        echo "[warn] ${tag}: attempt ${attempt}/3 failed (rc=${rc}); backing off"
        sleep $((attempt * 60))
    done
    if [ ${ok} -ne 1 ]; then
        echo "[FAIL] ${tag}: giving up after 3 attempts -- continuing with the rest"
        FAILED="${FAILED:-} ${tag}"
    fi
done

echo
echo "outputs under: ${OUTDIR}"
if [ -n "${FAILED:-}" ]; then
    echo "FAILED mass points:${FAILED}"
else
    echo "all requested mass points completed"
fi
