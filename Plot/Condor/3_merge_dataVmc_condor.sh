#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PLOT_DIR="${PROJECT_DIR}/Plot"
SCRIPTS_DIR="${PLOT_DIR}/scripts"
OUTPUT_DIR="${PLOT_DIR}/plots"
VARIABLES_DIR="${OUTPUT_DIR}/variables_dataVmc"
LOG_DIR="${OUTPUT_DIR}/logs_split"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RUN_DATAVMC_PLOTS="${RUN_DATAVMC_PLOTS:-1}"
RUN_OPTIMIZATION="${RUN_OPTIMIZATION:-1}"
RUN_MERGE_PLOTS="${RUN_MERGE_PLOTS:-1}"
CHECK_ROOT_KEYS="${CHECK_ROOT_KEYS:-1}"

sample_tags=(
    data
    dyll
    dyg
    sig_m1
    sig_m2
    sig_m3
    sig_m4
    sig_m5
    sig_m6
    sig_m7
    sig_m8
    sig_m9
    sig_m10
    sig_m15
    sig_m20
    sig_m25
    sig_m30
)

mkdir -p "$LOG_DIR" "$VARIABLES_DIR"
export PYTHONPATH="${PYTHONPATH:-}:${PLOT_DIR}/lib:${PROJECT_DIR}/HZaMVA/scripts"

cd "$PLOT_DIR"

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

root_has_keys() {
    local path="$1"

    if [[ "$CHECK_ROOT_KEYS" != "1" ]]; then
        return 0
    fi

    "$PYTHON_BIN" - "$path" <<'PY' >/dev/null 2>&1
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
}

merge_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local suffix target input sample_tag
    local inputs=()

    suffix="$(region_suffix "$region_key")"
    target="${VARIABLES_DIR}/ALP_plot_run3_${suffix}_${final_tag}.root"

    for sample_tag in "${sample_tags[@]}"; do
        input="${VARIABLES_DIR}/ALP_plot_run3_${suffix}_${final_tag}_part_${sample_tag}.root"
        if [[ ! -s "$input" ]]; then
            echo "[ERROR] Missing partial ROOT file: $input" >&2
            return 1
        fi
        if ! root_has_keys "$input"; then
            echo "[ERROR] Partial ROOT file has no keys or is unreadable: $input" >&2
            return 1
        fi
        inputs+=("$input")
    done

    echo "[hadd] $target"
    hadd -f "$target" "${inputs[@]}"
    if ! root_has_keys "$target"; then
        echo "[ERROR] Merged ROOT file has no keys or is unreadable: $target" >&2
        return 1
    fi
}

draw_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local log_file="${LOG_DIR}/draw_${final_tag}_${region_key}.log"
    local cmd=(
        "$PYTHON_BIN" "$SCRIPTS_DIR/2_plot_dataVmc.py"
        -y run3
        -m
        --ln
        --inputTag "$final_tag"
        --outputTag "$final_tag"
    )

    case "$region_key" in
        SR)  cmd+=(--region 1) ;;
        CR)  cmd+=(--region 2) ;;
        mva) cmd+=(-b) ;;
        *)
            echo "[ERROR] Unknown region key: $region_key" >&2
            return 1
            ;;
    esac

    {
        echo "[START] $(date '+%F %T') draw tag=${final_tag} region=${region_key}"
        printf "[CMD]"
        printf " %q" "${cmd[@]}"
        echo
    } > "$log_file"

    "${cmd[@]}" >> "$log_file" 2>&1
    echo "[DONE] $(date '+%F %T')" >> "$log_file"
}

if [[ "$RUN_MERGE_PLOTS" == "1" ]]; then
    merge_pids=()
    merge_labels=()

    for final_tag in nominal sideband_rwgt; do
        for region_key in SR CR mva; do
            echo "[submit hadd] tag=${final_tag} region=${region_key}"
            merge_plot_output "$region_key" "$final_tag" &
            merge_pids+=("$!")
            merge_labels+=("${final_tag}/${region_key}")
        done
    done

    merge_status=0
    for i in "${!merge_pids[@]}"; do
        if wait "${merge_pids[$i]}"; then
            echo "[done hadd] ${merge_labels[$i]}"
        else
            echo "[ERROR] hadd failed: ${merge_labels[$i]}" >&2
            merge_status=1
        fi
    done

    if [[ "$merge_status" != "0" ]]; then
        exit "$merge_status"
    fi
else
    echo "[Info] RUN_MERGE_PLOTS=$RUN_MERGE_PLOTS; skip hadd merge"
fi

if [[ "$RUN_DATAVMC_PLOTS" == "1" ]]; then
    for final_tag in nominal sideband_rwgt; do
        for region_key in SR CR mva; do
            draw_plot_output "$region_key" "$final_tag"
        done
    done
else
    echo "[Info] RUN_DATAVMC_PLOTS=$RUN_DATAVMC_PLOTS; skip 2_plot_dataVmc.py"
fi

if [[ "$RUN_OPTIMIZATION" == "1" ]]; then
    for final_tag in nominal sideband_rwgt; do
        "$PYTHON_BIN" "$SCRIPTS_DIR/ALP_Optimization.py" \
            -y run3 \
            -o "${OUTPUT_DIR}/optimize_run3UL_${final_tag}" \
            --region 1 \
            -p --sigVSscore -s --doOpt -c 1 \
            --inputTag "$final_tag"
    done
else
    echo "[Info] RUN_OPTIMIZATION=$RUN_OPTIMIZATION; skip ALP_Optimization.py"
fi
