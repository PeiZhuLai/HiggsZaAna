#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PLOT_DIR="${PROJECT_DIR}/Plot"
OUTPUT_DIR="${PLOT_DIR}/plots"
VARIABLES_DIR="${OUTPUT_DIR}/variables_dataVmc"

submit_file="${script_dir}/dataVmc.submit"
jobs_file="${script_dir}/dataVmc_jobs.txt"
resubmit_file="${script_dir}/dataVmc_resubmit.submit"
resubmit_jobs_file="${script_dir}/dataVmc_resubmit_jobs.txt"
missing_list="${script_dir}/dataVmc_missing_outputs.txt"

DRY_RUN="${DRY_RUN:-0}"
REMAKE_JOBS="${REMAKE_JOBS:-1}"

region_suffix() {
    local region_key="$1"
    case "$region_key" in
        SR)  echo "UL_SR" ;;
        CR)  echo "UL_CR" ;;
        mva) echo "UL_mva" ;;
        *)
            echo "[ERROR] Unknown region key: $region_key" >&2
            return 1
            ;;
    esac
}

expected_output() {
    local region_key="$1"
    local final_tag="$2"
    local sample_tag="$3"
    local suffix

    suffix="$(region_suffix "$region_key")"
    echo "${VARIABLES_DIR}/ALP_plot_run3_${suffix}_${final_tag}_part_${sample_tag}.root"
}

if [[ "$REMAKE_JOBS" == "1" || ! -s "$submit_file" || ! -s "$jobs_file" ]]; then
    if [[ "${SKIP_ENV_PACK:-0}" != "1" ]]; then
        "$script_dir/pack_conda_env_for_condor.sh"
    fi
    python3 "$script_dir/make_dataVmc_condor_jobs.py" --condor-dir "$script_dir"
fi

if [[ ! -s "$submit_file" ]]; then
    echo "[ERROR] Missing submit file: $submit_file" >&2
    exit 1
fi

if [[ ! -s "$jobs_file" ]]; then
    echo "[ERROR] Missing jobs file: $jobs_file" >&2
    exit 1
fi

mkdir -p "$VARIABLES_DIR"
: > "$resubmit_jobs_file"
: > "$missing_list"

total_jobs=0
missing_jobs=0

while read -r region_key final_tag sample_tag samples extra; do
    if [[ -z "${region_key:-}" || "${region_key:0:1}" == "#" ]]; then
        continue
    fi
    if [[ -n "${extra:-}" ]]; then
        echo "[ERROR] Unexpected extra fields in $jobs_file: $region_key $final_tag $sample_tag $samples $extra" >&2
        exit 1
    fi

    total_jobs=$((total_jobs + 1))
    output_path="$(expected_output "$region_key" "$final_tag" "$sample_tag")"

    if [[ ! -s "$output_path" ]]; then
        missing_jobs=$((missing_jobs + 1))
        printf "%s %s %s %s\n" "$region_key" "$final_tag" "$sample_tag" "$samples" >> "$resubmit_jobs_file"
        printf "%s %s %s %s -> %s\n" "$region_key" "$final_tag" "$sample_tag" "$samples" "$output_path" >> "$missing_list"
    fi
done < "$jobs_file"

echo "[Check] total jobs: $total_jobs"
echo "[Check] missing outputs: $missing_jobs"

if [[ "$missing_jobs" -eq 0 ]]; then
    rm -f "$resubmit_jobs_file" "$missing_list" "$resubmit_file"
    echo "[OK] All expected dataVmc partial ROOT files exist."
    exit 0
fi

echo "[Info] Missing output list: $missing_list"
echo "[Info] Resubmit jobs file: $resubmit_jobs_file"

awk -v jobs_file="$resubmit_jobs_file" '
    /^queue[[:space:]]+region_key,[[:space:]]*final_tag,[[:space:]]*sample_tag,[[:space:]]*samples[[:space:]]+from[[:space:]]+/ {
        print "queue region_key, final_tag, sample_tag, samples from " jobs_file
        next
    }
    { print }
' "$submit_file" > "$resubmit_file"

echo "[Info] Resubmit file: $resubmit_file"

if [[ "$DRY_RUN" == "1" ]]; then
    echo "[DryRun] Missing jobs that would be resubmitted:"
    sed -n '1,120p' "$resubmit_jobs_file"
    if [[ "$missing_jobs" -gt 120 ]]; then
        echo "[DryRun] ... truncated after 120 jobs"
    fi
    echo "[DryRun] Command: condor_submit $resubmit_file"
    exit 0
fi

condor_submit "$resubmit_file"
