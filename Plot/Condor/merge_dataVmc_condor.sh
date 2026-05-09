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
        inputs+=("$input")
    done

    echo "[hadd] $target"
    hadd -f "$target" "${inputs[@]}"
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

for final_tag in nominal sideband_rwgt; do
    for region_key in SR CR mva; do
        merge_plot_output "$region_key" "$final_tag"
    done
done

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
            --region 2 \
            -p --sigVSscore -s --doOpt -c 2 \
            --inputTag "$final_tag"
    done
else
    echo "[Info] RUN_OPTIMIZATION=$RUN_OPTIMIZATION; skip ALP_Optimization.py"
fi
