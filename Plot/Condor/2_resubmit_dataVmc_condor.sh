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
REMAKE_JOBS="${REMAKE_JOBS:-0}"
CHECK_ROOT_KEYS="${CHECK_ROOT_KEYS:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RESUBMIT_REQUEST_MEMORY="${RESUBMIT_REQUEST_MEMORY:-}"
RESUBMIT_REQUEST_DISK="${RESUBMIT_REQUEST_DISK:-}"
RESUBMIT_JOB_FLAVOUR="${RESUBMIT_JOB_FLAVOUR:-}"

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

output_problem() {
    local path="$1"

    if [[ ! -e "$path" ]]; then
        echo "missing"
        return 0
    fi

    if [[ ! -s "$path" ]]; then
        echo "empty"
        return 0
    fi

    if [[ "$CHECK_ROOT_KEYS" != "1" ]]; then
        return 1
    fi

    if "$PYTHON_BIN" - "$path" <<'PY' >/dev/null 2>&1; then
import sys

try:
    import ROOT
except Exception:
    sys.exit(2)

path = sys.argv[1]
ROOT.gROOT.SetBatch(True)
root_file = ROOT.TFile.Open(path, "READ")
if not root_file or root_file.IsZombie():
    sys.exit(1)

keys = root_file.GetListOfKeys()
sys.exit(0 if keys and keys.GetEntries() > 0 else 1)
PY
        return 1
    fi

    echo "no_root_keys_or_unreadable"
    return 0
}

if [[ "$REMAKE_JOBS" == "1" || ! -s "$submit_file" || ! -s "$jobs_file" ]]; then
    echo "[Info] Regenerating $submit_file and $jobs_file"
    if [[ "${SKIP_ENV_PACK:-0}" != "1" ]]; then
        "$script_dir/pack_conda_env_for_condor.sh"
    fi
    make_jobs_args=()
    if [[ -n "${DATA_VMC_MAKE_JOBS_ARGS:-}" ]]; then
        # shellcheck disable=SC2206
        make_jobs_args=(${DATA_VMC_MAKE_JOBS_ARGS})
    fi
    python3 "$script_dir/make_dataVmc_condor_jobs.py" --condor-dir "$script_dir" "${make_jobs_args[@]}"
else
    echo "[Info] Reusing existing submit template: $submit_file"
    echo "[Info] Reusing existing job list: $jobs_file"
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

    if problem="$(output_problem "$output_path")"; then
        missing_jobs=$((missing_jobs + 1))
        printf "%s %s %s %s\n" "$region_key" "$final_tag" "$sample_tag" "$samples" >> "$resubmit_jobs_file"
        printf "%s %s %s %s -> %s (%s)\n" "$region_key" "$final_tag" "$sample_tag" "$samples" "$output_path" "$problem" >> "$missing_list"
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

if [[ -n "$RESUBMIT_REQUEST_MEMORY" ]]; then
    awk -v value="$RESUBMIT_REQUEST_MEMORY" '
        /^request_memory[[:space:]]*=/ { print "request_memory = " value; next }
        { print }
    ' "$resubmit_file" > "${resubmit_file}.tmp"
    mv "${resubmit_file}.tmp" "$resubmit_file"
fi

if [[ -n "$RESUBMIT_REQUEST_DISK" ]]; then
    awk -v value="$RESUBMIT_REQUEST_DISK" '
        /^request_disk[[:space:]]*=/ { print "request_disk = " value; next }
        { print }
    ' "$resubmit_file" > "${resubmit_file}.tmp"
    mv "${resubmit_file}.tmp" "$resubmit_file"
fi

if [[ -n "$RESUBMIT_JOB_FLAVOUR" ]]; then
    awk -v value="$RESUBMIT_JOB_FLAVOUR" '
        /^\+JobFlavour[[:space:]]*=/ { print "+JobFlavour = \"" value "\""; next }
        { print }
    ' "$resubmit_file" > "${resubmit_file}.tmp"
    mv "${resubmit_file}.tmp" "$resubmit_file"
fi

echo "[Info] Resubmit file: $resubmit_file"
grep -E '^(request_memory|request_disk|\+JobFlavour)[[:space:]]*=' "$resubmit_file"

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
