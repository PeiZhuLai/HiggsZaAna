#!/bin/bash
# Submit MLNanoAOD production for all 4 DY background tags across all Run3 eras
# (read from HiggsDNA's zgamma_tutorial.json _MINI entries).
#
# Covers HiggsDNA tags:
#   DYGto2LG_10to100_MINI (2023preBPix / 2023postBPix / 2024)
#   DYGto2LG_10to50_MINI  (2022preEE / 2022postEE)
#   DYGto2LG_50to100_MINI (2022preEE / 2022postEE)
#   DYJetsToLL_MINI       (2022preEE / 2022postEE / 2023preBPix / 2023postBPix)
#   DYJetsTo2E_MINI       (2024) — replaces DYJetsToLL for 2024
#   DYJetsTo2Mu_MINI      (2024)
#   DYJetsTo2Tau_MINI     (2024)
#
# Each (tag, era) pair becomes one stage dir like Bkg_<tag>_<era>/ and one condor
# cluster. Outputs land at /eos/.../MLNanoAOD/Bkg_<tag>_<era>/.
#
# Usage:
#   submit_all_bkg.sh                       # ALL Run3 eras
#   submit_all_bkg.sh 2022preEE             # only that era
#   DRY_RUN=1 submit_all_bkg.sh             # print but don't submit
#   FILES_PER_JOB=8 submit_all_bkg.sh       # chunk 8 files per condor job
#
# Skip a (tag, era): set SKIP_DONE=1 (default) — if stage/<dir>/job.jdl exists
# it will be rebuilt but only re-submitted if EOS output dir is empty.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SAMPLE_JSON="${REPO_DIR}/HiggsDNA/metadata/samples/zgamma_tutorial.json"
SUBMITTER="${SCRIPT_DIR}/submit_mass_point.sh"

WANTED_ERA="${1:-all}"
DRY_RUN="${DRY_RUN:-0}"
FILES_PER_JOB="${FILES_PER_JOB:-8}"
SKIP_DONE="${SKIP_DONE:-1}"
export FILES_PER_JOB

BKG_TAGS=(
    DYGto2LG_10to100_MINI
    DYGto2LG_10to50_MINI
    DYGto2LG_50to100_MINI
    DYJetsToLL_MINI
    DYJetsTo2E_MINI
    DYJetsTo2Mu_MINI
    DYJetsTo2Tau_MINI
)

# Build (era, tag, das) tuples
mapfile -t entries < <(python3 - "${SAMPLE_JSON}" "${WANTED_ERA}" "${BKG_TAGS[@]}" <<'PY'
import json, sys
sample_json, wanted_era, *tags = sys.argv[1:]
with open(sample_json) as f:
    d = json.load(f)
RUN3_ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
eras = RUN3_ERAS if wanted_era == "all" else [wanted_era]
for tag in tags:
    if tag not in d:
        continue
    files = d[tag].get("files", {})
    for era in eras:
        paths = files.get(era, [])
        if not paths:
            continue
        # Use only the first path (these are MINIAOD parent datasets, normally singletons)
        path = paths[0].strip() if isinstance(paths[0], str) else ""
        if not path or not path.startswith("/"):
            continue
        # tag-without-_MINI for shorter stage dirs
        short = tag.replace("_MINI", "")
        out_tag = f"Bkg_{short}_{era}"
        print(f"{era}|{out_tag}|{path}")
PY
)

echo "Found ${#entries[@]} (tag, era) bkg datasets to consider"
echo

total_submitted=0
total_skipped=0
for entry in "${entries[@]}"; do
    IFS='|' read -r era out_tag das <<< "${entry}"
    echo "--- ${out_tag} ---"
    echo "  era=${era}"
    echo "  das=${das}"

    # skip if already finished (EOS dir non-empty)
    eos_dir="/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/${out_tag}"
    n_done=$(ls "${eos_dir}" 2>/dev/null | wc -l)
    if [ "${SKIP_DONE}" = "1" ] && [ "${n_done}" -gt 0 ]; then
        echo "  SKIP: EOS already has ${n_done} files at ${eos_dir}"
        total_skipped=$((total_skipped + 1))
        continue
    fi

    if [ "${DRY_RUN}" = "1" ]; then
        continue
    fi

    # prep stage (writes files.txt, args.txt, JDL)
    bash "${SUBMITTER}" "${out_tag}" "${das}" mc "${era}" 2>&1 | tail -3
    jdl="${REPO_DIR}/RegressMergedPhoton/condor/stage/${out_tag}/job.jdl"
    result=$(condor_submit "${jdl}" 2>&1 | tail -1)
    echo "  ${result}"
    total_submitted=$((total_submitted + 1))
done

echo
echo "Submitted ${total_submitted} clusters, skipped ${total_skipped} (already done)."
