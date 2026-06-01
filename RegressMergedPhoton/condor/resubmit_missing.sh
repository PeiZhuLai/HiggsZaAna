#!/bin/bash
# Re-submit only the still-missing output groups for every stage/<tag>/ that has
# an args.txt. Each output ROOT name is checked against the EOS output dir; if
# the file already exists, the corresponding args line is skipped.
#
# Usage:
#   resubmit_missing.sh                  # all stage dirs
#   resubmit_missing.sh 'mA_M*'          # only matching tag (glob, quoted)
#   DRY_RUN=1 resubmit_missing.sh ...    # print what would be submitted
#
# Produces, for each tag with anything missing:
#   stage/<tag>/args_missing.txt
#   stage/<tag>/job_missing.jdl
# and condor_submits job_missing.jdl.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
STAGE_ROOT="${REPO_DIR}/RegressMergedPhoton/condor/stage"

FILTER="${1:-*}"
DRY_RUN="${DRY_RUN:-0}"

echo "Stage root  : ${STAGE_ROOT}"
echo "Filter      : ${FILTER}"
echo "DRY_RUN     : ${DRY_RUN}"
echo

total_clusters=0
total_jobs=0
total_skipped=0

shopt -s nullglob
for stage in "${STAGE_ROOT}"/${FILTER}/; do
    tag=$(basename "${stage}")
    args="${stage}/args.txt"
    jdl="${stage}/job.jdl"
    [ -f "${args}" ] || { echo "[SKIP] ${tag}: no args.txt"; continue; }
    [ -f "${jdl}"  ] || { echo "[SKIP] ${tag}: no job.jdl";  continue; }

    args_missing="${stage}/args_missing.txt"
    awk '{
        # cols: input_urls out_name out_dir
        cmd = "test -f " $3 "/" $2;
        if (system(cmd) != 0) print $0;
    }' "${args}" > "${args_missing}"

    n_total=$(wc -l < "${args}")
    n_miss=$(wc -l < "${args_missing}")
    n_done=$((n_total - n_miss))

    if [ "${n_miss}" -eq 0 ]; then
        echo "[DONE] ${tag}: ${n_total}/${n_total} files already on EOS"
        total_skipped=$((total_skipped + 1))
        rm -f "${args_missing}"
        continue
    fi

    jdl_missing="${stage}/job_missing.jdl"
    sed "s|${args}|${args_missing}|" "${jdl}" > "${jdl_missing}"

    echo "[SUBMIT] ${tag}: ${n_miss}/${n_total} missing (${n_done} done)"
    if [ "${DRY_RUN}" = "1" ]; then
        continue
    fi

    result=$(condor_submit "${jdl_missing}" 2>&1 | tail -1)
    echo "    ${result}"
    total_clusters=$((total_clusters + 1))
    total_jobs=$((total_jobs + n_miss))
done

echo
echo "Summary: submitted ${total_clusters} clusters (${total_jobs} jobs); ${total_skipped} tags already complete."
