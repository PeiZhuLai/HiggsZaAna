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
SETUP_CONDA_ENV="${SETUP_CONDA_ENV:-auto}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
LOCALIZE_CONDA_ENV="${LOCALIZE_CONDA_ENV:-auto}"
LOCAL_CONDA_DIR="${LOCAL_CONDA_DIR:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/higgs-alp-ana-conda}"
CONDA_TARBALL="${CONDA_TARBALL:-${PROJECT_DIR}/Plot/Condor/env_cache/higgs-alp-ana.tar.gz}"

activate_conda_env_if_needed() {
    local should_activate=0

    case "$SETUP_CONDA_ENV" in
        1|true|TRUE|yes|YES) should_activate=1 ;;
        0|false|FALSE|no|NO) should_activate=0 ;;
        auto|AUTO)
            if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV_NAME" ]]; then
                should_activate=1
            fi
            ;;
        *)
            echo "[ERROR] SETUP_CONDA_ENV must be auto, 1, or 0; got '$SETUP_CONDA_ENV'" >&2
            exit 2
            ;;
    esac

    if [[ "$should_activate" -ne 1 ]]; then
        return 0
    fi

    echo "[ENV] Activate conda env: $CONDA_ENV_NAME"
    if [[ ! -r "$ANACONDA_SETUP" ]]; then
        echo "[ERROR] Cannot read ANACONDA_SETUP: $ANACONDA_SETUP" >&2
        exit 2
    fi

    echo "[ENV] Source anaconda setup: $ANACONDA_SETUP"
    set +u
    source "$ANACONDA_SETUP"
    conda activate "$CONDA_ENV_NAME"
    set -u

    if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/python3" ]]; then
        echo "[ERROR] Failed to activate conda env: $CONDA_ENV_NAME" >&2
        exit 2
    fi

    PYTHON_BIN="${CONDA_PREFIX}/bin/python3"
    echo "[ENV] active CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
    echo "[ENV] active python=${PYTHON_BIN}"
}

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

    if [[ ! -s "$CONDA_TARBALL" ]]; then
        echo "[ERROR] Conda localization requested, but tarball is missing: $CONDA_TARBALL" >&2
        echo "[ERROR] Run Plot/Condor/pack_conda_env_for_condor.sh once before condor_submit." >&2
        exit 2
    fi

    echo "[ENV] Extract conda tarball from $CONDA_TARBALL to $LOCAL_CONDA_DIR"
    rm -rf "$LOCAL_CONDA_DIR"
    mkdir -p "$LOCAL_CONDA_DIR"
    tar -xzf "$CONDA_TARBALL" -C "$LOCAL_CONDA_DIR"

    if [[ -x "$LOCAL_CONDA_DIR/bin/conda-unpack" ]]; then
        echo "[ENV] Run conda-unpack"
        "$LOCAL_CONDA_DIR/bin/conda-unpack"
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
    echo "[ENV] SETUP_CONDA_ENV=${SETUP_CONDA_ENV}"
    echo "[ENV] CONDA_ENV_NAME=${CONDA_ENV_NAME}"
    echo "[ENV] ANACONDA_SETUP=${ANACONDA_SETUP}"
    echo "[ENV] LOCALIZE_CONDA_ENV=${LOCALIZE_CONDA_ENV}"
    echo "[ENV] CONDA_TARBALL=${CONDA_TARBALL}"
} > "$log_file"

activate_conda_env_if_needed >> "$log_file" 2>&1
localize_conda_env_if_needed >> "$log_file" 2>&1
"$PYTHON_BIN" -c 'import numpy; import ROOT; import xgboost; print("[ENV] python imports OK")' >> "$log_file" 2>&1

cmd=(
    "$PYTHON_BIN" "$SCRIPTS_DIR/1_prepare_dataVmc.py"
    -y run3
    -m
    --ln
    --histOnly
    --skipSystematics
    --optimizeBranches
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
