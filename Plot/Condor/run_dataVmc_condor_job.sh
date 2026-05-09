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
LOCALIZE_CONDA_ENV="${LOCALIZE_CONDA_ENV:-auto}"
LOCAL_CONDA_DIR="${LOCAL_CONDA_DIR:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/higgs-alp-ana-conda}"

localize_conda_env_if_needed() {
    local source_env="${CONDA_PREFIX:-}"
    local should_localize=0

    case "$LOCALIZE_CONDA_ENV" in
        1|true|TRUE|yes|YES) should_localize=1 ;;
        0|false|FALSE|no|NO) should_localize=0 ;;
        auto|AUTO)
            if [[ -n "$source_env" && "$source_env" == /eos/* ]]; then
                should_localize=1
            fi
            ;;
        *)
            echo "[ERROR] LOCALIZE_CONDA_ENV must be auto, 1, or 0; got '$LOCALIZE_CONDA_ENV'" >&2
            exit 2
            ;;
    esac

    if [[ "$should_localize" -ne 1 ]]; then
        return 0
    fi
    if [[ -z "$source_env" || ! -x "$source_env/bin/python3" ]]; then
        echo "[WARN] LOCALIZE_CONDA_ENV requested, but CONDA_PREFIX is not a usable conda env: ${source_env:-<unset>}" >&2
        return 0
    fi

    echo "[ENV] Copy conda env from $source_env to $LOCAL_CONDA_DIR"
    mkdir -p "$(dirname "$LOCAL_CONDA_DIR")"
    if command -v rsync >/dev/null 2>&1; then
        rsync -a --delete \
            --exclude '__pycache__/' \
            --exclude '*.pyc' \
            "$source_env/" "$LOCAL_CONDA_DIR/"
    else
        rm -rf "$LOCAL_CONDA_DIR"
        cp -a "$source_env" "$LOCAL_CONDA_DIR"
    fi

    export CONDA_PREFIX="$LOCAL_CONDA_DIR"
    export PATH="$LOCAL_CONDA_DIR/bin:$PATH"
    export LD_LIBRARY_PATH="$LOCAL_CONDA_DIR/lib:${LD_LIBRARY_PATH:-}"
    PYTHON_BIN="$LOCAL_CONDA_DIR/bin/python3"
}

mkdir -p "$LOG_DIR" "$VARIABLES_DIR"
export PYTHONPATH="${PYTHONPATH:-}:${PLOT_DIR}/lib:${PROJECT_DIR}/HZaMVA/scripts"

partial_tag="${final_tag}_part_${sample_tag}"
log_file="${LOG_DIR}/${final_tag}_${region_key}_${sample_tag}.log"

cd "$PLOT_DIR"

{
    echo "[START] $(date '+%F %T')"
    echo "[HOST] $(hostname)"
    echo "[PWD] $(pwd)"
    echo "[JOB] tag=${final_tag} region=${region_key} sample_tag=${sample_tag} samples=${samples}"
    echo "[ENV] initial PYTHON_BIN=${PYTHON_BIN}"
    echo "[ENV] initial CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
    echo "[ENV] LOCALIZE_CONDA_ENV=${LOCALIZE_CONDA_ENV}"
} > "$log_file"

localize_conda_env_if_needed >> "$log_file" 2>&1
"$PYTHON_BIN" -c 'import numpy; import ROOT; import xgboost; print("[ENV] python imports OK")' >> "$log_file" 2>&1

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
    printf "[CMD]"
    printf " %q" "${cmd[@]}"
    echo
} >> "$log_file"

"${cmd[@]}" >> "$log_file" 2>&1
echo "[DONE] $(date '+%F %T')" >> "$log_file"
