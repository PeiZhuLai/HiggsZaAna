#!/bin/bash
# Like submit_all_data.sh, but SKIP datasets whose MLNanoAOD output is already
# complete (output .root count >= #files in dataset). Useful for:
#   - submitting the remaining 31 datasets after a pilot dataset is done
#   - resubmitting only the incomplete datasets after held/failed jobs
#
# The worker writes ONE output per input MiniAOD file (named by UUID), so a
# dataset is "complete" when (# outputs in OUT_EOS) >= (# files in dataset).
#
# Usage:
#   submit_remaining_data.sh 2024            # 2024 only (default era arg)
#   DRY_RUN=1 submit_remaining_data.sh 2024  # show plan, don't submit
#   FILES_PER_JOB=8 submit_remaining_data.sh 2024
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SAMPLE_JSON="${REPO_DIR}/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json"
[ -f "${SAMPLE_JSON}" ] || SAMPLE_JSON="$(cd "${REPO_DIR}/.." && pwd)/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json"
SUBMITTER="${SCRIPT_DIR}/submit_mass_point.sh"
OUT_BASE="/eos/home-p/pelai/HZa/MLNanoAOD"

WANTED_ERA="${1:-2024}"
DRY_RUN="${DRY_RUN:-0}"
FILES_PER_JOB="${FILES_PER_JOB:-8}"
export FILES_PER_JOB

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
        m = re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+", path)
        if not m:
            continue
        primary, _yr, run, rest = m.group(1), m.group(2), m.group(3), m.group(4)
        vers = "".join(re.findall(r"[-_]v\d+", rest)).replace("-", "").replace("_", "")
        tag = f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}"
        print(f"{era}|{tag}|{path}")
PY
)

echo "Found ${#entries[@]} data datasets for era=${WANTED_ERA}"
echo "FILES_PER_JOB=${FILES_PER_JOB}  OUT_BASE=${OUT_BASE}"
echo

source /cvmfs/cms.cern.ch/cmsset_default.sh > /dev/null 2>&1
submitted=0; skipped=0
for entry in "${entries[@]}"; do
    IFS='|' read -r era tag das <<< "${entry}"

    # nfiles: prefer cached files.txt from a prior prep, else query DAS
    ftxt="${SCRIPT_DIR}/stage/${tag}/files.txt"
    if [ -s "${ftxt}" ]; then
        nfiles=$(wc -l < "${ftxt}")
    else
        nfiles=$(dasgoclient -query "file dataset=${das}" 2>/dev/null | wc -l)
    fi
    nout=$(ls "${OUT_BASE}/${tag}"/*.root 2>/dev/null | wc -l)

    if [ "${nfiles}" -gt 0 ] && [ "${nout}" -ge "${nfiles}" ]; then
        echo "SKIP  ${tag}  (complete: ${nout}/${nfiles} outputs)"
        skipped=$((skipped+1))
        continue
    fi

    echo "--- ${tag}  (have ${nout}/${nfiles}) ---"
    if [ "${DRY_RUN}" = "1" ]; then continue; fi

    bash "${SUBMITTER}" "${tag}" "${das}" data "${era}" 2>&1 | tail -3
    jdl="${SCRIPT_DIR}/stage/${tag}/job.jdl"
    result=$(condor_submit "${jdl}" 2>&1 | tail -1)
    echo "  ${result}"
    submitted=$((submitted+1))
done
echo
echo "Done. submitted=${submitted} datasets, skipped=${skipped} (already complete)."
