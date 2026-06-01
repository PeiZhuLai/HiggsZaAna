#!/bin/bash
# Submit MLNanoAOD production for all Run3 data MiniAOD datasets defined in
# zgamma_tutorial_sample_manager_full.json under the "Data_MINI" key.
#
# Per dataset, a tag like "EGamma_2024_C" is built and submitted as one cluster.
# Outputs land at /eos/.../MLNanoAOD/Data_<era>_<primary>_<run>/
#
# Usage:
#   submit_all_data.sh                # submit ALL Run3 eras
#   submit_all_data.sh 2024           # only specific era
#   DRY_RUN=1 submit_all_data.sh      # print but don't submit
#
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SAMPLE_JSON="${REPO_DIR}/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json"
SUBMITTER="${SCRIPT_DIR}/submit_mass_point.sh"

WANTED_ERA="${1:-all}"
DRY_RUN="${DRY_RUN:-0}"
FILES_PER_JOB="${FILES_PER_JOB:-8}"
export FILES_PER_JOB

# Use python to list (era, dataset) pairs from the JSON
mapfile -t entries < <(python3 - "${SAMPLE_JSON}" "${WANTED_ERA}" <<'PY'
import json, sys, re
with open(sys.argv[1]) as f:
    d = json.load(f)
wanted = sys.argv[2]
eras = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
if wanted != "all":
    eras = [e for e in eras if e == wanted]
for era in eras:
    for path in d["Data_MINI"]["files"].get(era, []):
        # /EGamma0/Run2024G-MINIv6NANOv15-v2/MINIAOD -> primary=EGamma0, run=G, ver=v2
        # /EGamma1/Run2023C-22Sep2023_v3-v1/MINIAOD  -> primary=EGamma1, run=C, ver=v3v1
        m = re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+", path)
        if not m:
            continue
        primary, _yr, run, rest = m.group(1), m.group(2), m.group(3), m.group(4)
        # extract all _vN or -vN tokens from rest, join → ver string
        vers = "".join(re.findall(r"[-_]v\d+", rest)).replace("-", "").replace("_", "")
        tag = f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}"
        print(f"{era}|{tag}|{path}")
PY
)

echo "Found ${#entries[@]} data datasets to submit"
echo

total_clusters=0
for entry in "${entries[@]}"; do
    IFS='|' read -r era tag das <<< "${entry}"
    echo "--- ${tag} ---"
    echo "  era=${era}"
    echo "  das=${das}"
    if [ "${DRY_RUN}" = "1" ]; then
        continue
    fi
    # prep stage (writes files.txt, args.txt, JDL)
    bash "${SUBMITTER}" "${tag}" "${das}" data "${era}" 2>&1 | tail -5
    jdl="${REPO_DIR}/RegressMergedPhoton/condor/stage/${tag}/job.jdl"
    result=$(condor_submit "${jdl}" 2>&1 | tail -1)
    echo "  ${result}"
    total_clusters=$((total_clusters + 1))
done

echo
echo "Submitted ${total_clusters} clusters."
