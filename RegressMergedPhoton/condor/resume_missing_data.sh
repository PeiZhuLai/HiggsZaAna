#!/bin/bash
# Submit ONLY the missing MiniAOD files per dataset (output-presence is the ledger).
#
# The worker writes one output "<tag>_<UUID>.root" per input file. A file is
# "done" when that output exists in OUT_BASE/<tag>. This script, per dataset:
#   1. gets the full file list (cached stage/<tag>/files.txt, else dasgoclient),
#   2. computes which UUIDs have no output yet (= missing),
#   3. if none missing -> SKIP (complete),
#   4. else builds fresh groups of ONLY the missing URLs and condor_submits them.
#
# This subsumes a first-pass submission (fresh dataset = all files missing) AND
# straggler recovery (a few transient xrootd failures) without reprocessing
# already-done files.
#
# Usage:
#   resume_missing_data.sh 2024            # 2024 only (default)
#   DRY_RUN=1 resume_missing_data.sh 2024  # show counts, don't submit
#   FILES_PER_JOB=8 resume_missing_data.sh 2024
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SAMPLE_JSON="${REPO_DIR}/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json"
[ -f "${SAMPLE_JSON}" ] || SAMPLE_JSON="$(cd "${REPO_DIR}/.." && pwd)/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json"

WORKER="${SCRIPT_DIR}/run_one_file.sh"
PROXY="${X509_USER_PROXY:-/afs/cern.ch/user/p/pelai/.x509up_pelai}"
OUT_BASE="/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD"
WANTED_ERA="${1:-2024}"
DRY_RUN="${DRY_RUN:-0}"
FILES_PER_JOB="${FILES_PER_JOB:-8}"
MODE="data"

mapfile -t entries < <(python3 - "${SAMPLE_JSON}" "${WANTED_ERA}" <<'PY'
import json, sys, re
with open(sys.argv[1]) as f:
    d = json.load(f)
wanted = sys.argv[2]
eras = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
if wanted != "all":
    eras = [e for e in eras if e == wanted]
for era in eras:
    for path in d["Data_MINI"]["files"].get(era, []):
        m = re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+", path)
        if not m: continue
        primary, _yr, run, rest = m.group(1), m.group(2), m.group(3), m.group(4)
        vers = "".join(re.findall(r"[-_]v\d+", rest)).replace("-","").replace("_","")
        tag = f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}"
        print(f"{era}|{tag}|{path}")
PY
)

echo "era=${WANTED_ERA}  datasets=${#entries[@]}  FILES_PER_JOB=${FILES_PER_JOB}  OUT_BASE=${OUT_BASE}"
echo
source /cvmfs/cms.cern.ch/cmsset_default.sh > /dev/null 2>&1

total_jobs=0; total_missing=0; n_submit=0; n_skip=0
for entry in "${entries[@]}"; do
    IFS='|' read -r era tag das <<< "${entry}"
    STAGE="${SCRIPT_DIR}/stage/${tag}"
    OUT_EOS="${OUT_BASE}/${tag}"
    mkdir -p "${STAGE}"

    # full file list (cache)
    FILES_TXT="${STAGE}/files.txt"
    if [ ! -s "${FILES_TXT}" ]; then
        dasgoclient -query "file dataset=${das}" > "${FILES_TXT}"
    fi
    nfiles=$(wc -l < "${FILES_TXT}")

    # missing = input files whose UUID has no output yet. Match by UUID (the
    # final _-delimited token, since UUIDs have no underscores) so BOTH naming
    # schemes count as done: "<tag>_<uuid>.root" and "<tag>_resume_<uuid>.root".
    DONE_TXT="${STAGE}/done_uuids.txt"
    ls "${OUT_EOS}"/*.root 2>/dev/null | sed -E 's#.*_([0-9a-fA-F-]+)\.root$#\1#' | sort -u > "${DONE_TXT}"
    MISS_TXT="${STAGE}/missing.txt"
    awk 'NR==FNR{done[$0]=1;next}{n=$0; sub(/.*\//,"",n); sub(/\.root$/,"",n); if(!(n in done)) print "root://cms-xrd-global.cern.ch/"$0}' \
        "${DONE_TXT}" "${FILES_TXT}" > "${MISS_TXT}"
    nmiss=$(wc -l < "${MISS_TXT}")

    if [ "${nmiss}" -eq 0 ]; then
        echo "SKIP  ${tag}  (complete ${nfiles}/${nfiles})"
        n_skip=$((n_skip+1)); continue
    fi

    njobs=$(( (nmiss + FILES_PER_JOB - 1) / FILES_PER_JOB ))
    echo "RESUB ${tag}  missing ${nmiss}/${nfiles}  -> ${njobs} jobs"
    total_missing=$((total_missing+nmiss)); total_jobs=$((total_jobs+njobs))
    [ "${DRY_RUN}" = "1" ] && continue

    # build resume groups + args
    RDIR="${STAGE}/resume"
    rm -rf "${RDIR}"; mkdir -p "${RDIR}/inputs" "${RDIR}/logs"
    ARGS="${RDIR}/args.txt"; : > "${ARGS}"
    i=0; gi=0; gf=""
    while IFS= read -r url; do
        i=$((i+1))
        if [ -z "${gf}" ]; then gi=$((gi+1)); gf="${RDIR}/inputs/group_$(printf '%04d' ${gi}).txt"; : > "${gf}"; fi
        echo "${url}" >> "${gf}"
        if [ $((i % FILES_PER_JOB)) -eq 0 ]; then
            echo "${gf} ${tag}_resume_$(printf '%04d' ${gi}).root ${OUT_EOS}" >> "${ARGS}"; gf=""
        fi
    done < "${MISS_TXT}"
    [ -n "${gf}" ] && echo "${gf} ${tag}_resume_$(printf '%04d' ${gi}).root ${OUT_EOS}" >> "${ARGS}"

    JDL="${RDIR}/job.jdl"
    cat > "${JDL}" <<JDLEOF
Universe                = vanilla
Executable              = ${WORKER}
Arguments               = \$(input_url) \$(out_name) \$(out_dir) ${MODE} ${era}
should_transfer_files   = NO
use_x509userproxy       = true
x509userproxy           = ${PROXY}
JobBatchName            = "MLPhoton_resume_${tag}"
Log                     = ${RDIR}/logs/all_jobs.log
Output                  = ${RDIR}/logs/job_\$(Process).out
Error                   = ${RDIR}/logs/job_\$(Process).err
+JobFlavour             = "workday"
RequestCpus             = 1
RequestMemory           = 4000
RequestDisk             = 4000000
max_materialize         = 40
Queue input_url out_name out_dir from ${ARGS}
JDLEOF
    res=$(condor_submit "${JDL}" 2>&1 | tail -1)
    echo "  ${res}"
    n_submit=$((n_submit+1))
done
echo
echo "Summary: submitted ${n_submit} datasets (${total_jobs} jobs, ${total_missing} missing files); skipped ${n_skip} complete."
