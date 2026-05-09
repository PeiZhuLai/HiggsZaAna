#!/bin/bash
set -euo pipefail

timer_start=$(date +%s)
finish() {
    local status=$?
    local timer_end elapsed hours minutes seconds
    timer_end=$(date +%s)
    elapsed=$((timer_end - timer_start))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    printf "[RUNTIME] %02d:%02d:%02d\n" "$hours" "$minutes" "$seconds"
    exit "$status"
}
trap finish EXIT

if [[ "$#" -ne 4 ]]; then
    echo "[ERROR] Usage: $0 <region_key> <final_tag> <sample_tag> <samples>" >&2
    exit 2
fi

region_key="$1"
final_tag="$2"
sample_tag="$3"
samples="$4"

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PLOT_DIR="${PROJECT_DIR}/Plot"
SCRIPTS_DIR="${PLOT_DIR}/scripts"
OUTPUT_DIR="${PLOT_DIR}/plots"
VARIABLES_DIR="${OUTPUT_DIR}/variables_dataVmc"
LOG_DIR="${OUTPUT_DIR}/logs_split"
SIDEBAND_REWEIGHT_JSON="${SIDEBAND_REWEIGHT_JSON:-${PROJECT_DIR}/HZaMVA/reweights/sideband_run3_iterative.json}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MAX_EVENTS="${MAX_EVENTS:-}"

mkdir -p "$LOG_DIR" "$VARIABLES_DIR"
export PYTHONPATH="${PYTHONPATH:-}:${PLOT_DIR}/lib:${PROJECT_DIR}/HZaMVA/scripts"

cd "$PLOT_DIR"

partial_tag="${final_tag}_part_${sample_tag}"
log_file="${LOG_DIR}/${final_tag}_${region_key}_${sample_tag}.log"
cmd=(
    "$PYTHON_BIN" "$SCRIPTS_DIR/1_prepare_dataVmc.py"
    -y run3
    -m
    --ln
    --histOnly
    --samples "$samples"
    --outputTag "$partial_tag"
)

case "$region_key" in
    SR)  cmd+=(--region 1) ;;
    CR)  cmd+=(--region 2) ;;
    mva) cmd+=(-b) ;;
    *)
        echo "[ERROR] Unknown region key: $region_key" >&2
        exit 2
        ;;
esac

if [[ "$final_tag" == "sideband_rwgt" ]]; then
    cmd+=(--useSidebandReweight --sidebandReweightJson "$SIDEBAND_REWEIGHT_JSON")
fi

if [[ -n "$MAX_EVENTS" ]]; then
    cmd+=(--maxEvents "$MAX_EVENTS")
fi

{
    echo "[START] $(date '+%F %T')"
    echo "[HOST] $(hostname)"
    echo "[PWD] $(pwd)"
    echo "[JOB] tag=${final_tag} region=${region_key} sample_tag=${sample_tag} samples=${samples}"
    printf "[CMD]"
    printf " %q" "${cmd[@]}"
    echo
} > "$log_file"

"${cmd[@]}" >> "$log_file" 2>&1
echo "[DONE] $(date '+%F %T')" >> "$log_file"
