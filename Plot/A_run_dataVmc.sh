#!/bin/bash
timer_start=$(date +%s)
echo "==============TIMER STARTED: $(date '+%F %T')=============="
print_runtime() {
    local timer_end elapsed hours minutes seconds
    timer_end=$(date +%s)
    elapsed=$((timer_end - timer_start))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    printf "==============RUNTIME: %02d:%02d:%02d==============\n" "$hours" "$minutes" "$seconds"
}
trap print_runtime EXIT

pwd

# 0: Full Region, 1: Signal Region, 2: Contral Region
#-----------------------------Parallel Run-----------------------------------
set -euo pipefail

scriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts'
outputDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots'
variablesDir="${outputDir}/variables_dataVmc"
sidebandReweightJson='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json'
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib"

MAX_PLOT_JOBS="${MAX_PLOT_JOBS:-6}"
RUN_PREPARE_DATAVMC="${RUN_PREPARE_DATAVMC:-1}"
RUN_MERGE_PLOTS="${RUN_MERGE_PLOTS:-1}"
RUN_OPTIMIZATION="${RUN_OPTIMIZATION:-1}"
RUN_DATAVMC_PLOTS="${RUN_DATAVMC_PLOTS:-1}"
MAX_EVENTS="${MAX_EVENTS:-}"
logDir="${outputDir}/logs_split"
mkdir -p "$logDir" "$variablesDir"

sampleGroups=(
    "Data"
    "DYJetsToLL"
    "DYGto2LG"
    "M1,M2,M3,M4,M5"
    "M6,M7,M8,M9,M10"
    "M15,M20,M25,M30"
)

sampleTags=(
    "data"
    "dyll"
    "dyg"
    "sig_m1_m5"
    "sig_m6_m10"
    "sig_m15_m30"
)

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

run_plot_task() {
    local region_key="$1"
    local final_tag="$2"
    local sample_tag="$3"
    local samples="$4"
    local partial_tag="${final_tag}_part_${sample_tag}"
    local log_file="${logDir}/${final_tag}_${region_key}_${sample_tag}.log"
    local cmd=(python3 "$scriptsDir/1_prepare_dataVmc.py" -y run3 -m --ln --histOnly --samples "$samples" --outputTag "$partial_tag")

    case "$region_key" in
        SR)  cmd+=(--region 1) ;;
        CR)  cmd+=(--region 2) ;;
        mva) cmd+=(-b) ;;
        *)
            echo "[ERROR] Unknown region key: $region_key" >&2
            return 1
            ;;
    esac

    if [[ "$final_tag" == "sideband_rwgt" ]]; then
        cmd+=(--useSidebandReweight --sidebandReweightJson "$sidebandReweightJson")
    fi
    if [[ -n "$MAX_EVENTS" ]]; then
        cmd+=(--maxEvents "$MAX_EVENTS")
    fi

    {
        echo "[START] $(date '+%F %T') tag=${final_tag} region=${region_key} samples=${samples}"
        printf '[CMD]'
        printf ' %q' "${cmd[@]}"
        echo
    } > "$log_file"

    "${cmd[@]}" >> "$log_file" 2>&1
    echo "[DONE] $(date '+%F %T')" >> "$log_file"
}

if [[ "$RUN_PREPARE_DATAVMC" == "1" ]]; then
    activePids=()
    activeLabels=()
    plot_failed=0

    reap_finished_jobs() {
        local running_pids
        running_pids=" $(jobs -pr | tr '\n' ' ') "
        local kept_pids=()
        local kept_labels=()
        local idx pid label

        for idx in "${!activePids[@]}"; do
            pid="${activePids[$idx]}"
            label="${activeLabels[$idx]}"
            if [[ "$running_pids" == *" ${pid} "* ]]; then
                kept_pids+=("$pid")
                kept_labels+=("$label")
            else
                if ! wait "$pid"; then
                    echo "[ERROR] Plot job failed: $label. Check logs in: $logDir" >&2
                    plot_failed=1
                fi
            fi
        done

        activePids=("${kept_pids[@]}")
        activeLabels=("${kept_labels[@]}")
    }

    wait_for_plot_slot() {
        while true; do
            reap_finished_jobs
            if [[ "${#activePids[@]}" -lt "$MAX_PLOT_JOBS" ]]; then
                return 0
            fi
            sleep 10
        done
    }

    launch_plot_task() {
        wait_for_plot_slot
        run_plot_task "$@" &
        activePids+=("$!")
        activeLabels+=("tag=$2 region=$1 samples=$4")
    }

    for finalTag in nominal sideband_rwgt; do
        for regionKey in SR CR mva; do
            for i in "${!sampleGroups[@]}"; do
                launch_plot_task "$regionKey" "$finalTag" "${sampleTags[$i]}" "${sampleGroups[$i]}"
            done
        done
    done

    while [[ "${#activePids[@]}" -gt 0 ]]; do
        reap_finished_jobs
        if [[ "${#activePids[@]}" -gt 0 ]]; then
            sleep 10
        fi
    done

    if [[ "$plot_failed" -ne 0 ]]; then
        echo "[ERROR] At least one plot job failed. Check logs in: $logDir" >&2
        exit 1
    fi
else
    echo "[Info] RUN_PREPARE_DATAVMC=$RUN_PREPARE_DATAVMC; skip 1_prepare_dataVmc.py"
fi

merge_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local suffix
    suffix="$(region_suffix "$region_key")"
    local target="${variablesDir}/ALP_plot_run3_${suffix}_${final_tag}.root"
    local inputs=()

    for sample_tag in "${sampleTags[@]}"; do
        local input="${variablesDir}/ALP_plot_run3_${suffix}_${final_tag}_part_${sample_tag}.root"
        if [[ ! -s "$input" ]]; then
            echo "[ERROR] Missing partial ROOT file: $input" >&2
            return 1
        fi
        inputs+=("$input")
    done

    echo "[hadd] $target"
    hadd -f "$target" "${inputs[@]}"
}

if [[ "$RUN_MERGE_PLOTS" == "1" ]]; then
    for finalTag in nominal sideband_rwgt; do
        for regionKey in SR CR mva; do
            merge_plot_output "$regionKey" "$finalTag"
        done
    done
else
    echo "[Info] RUN_MERGE_PLOTS=$RUN_MERGE_PLOTS; skip hadd merge"
fi

draw_plot_output() {
    local region_key="$1"
    local final_tag="$2"
    local log_file="${logDir}/draw_${final_tag}_${region_key}.log"
    local cmd=(python3 "$scriptsDir/2_plot_dataVmc.py" -y run3 -m --ln --inputTag "$final_tag" --outputTag "$final_tag")

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
        printf '[CMD]'
        printf ' %q' "${cmd[@]}"
        echo
    } > "$log_file"

    "${cmd[@]}" >> "$log_file" 2>&1
    echo "[DONE] $(date '+%F %T')" >> "$log_file"
}

if [[ "$RUN_DATAVMC_PLOTS" == "1" ]]; then
    for finalTag in nominal sideband_rwgt; do
        for regionKey in SR CR mva; do
            draw_plot_output "$regionKey" "$finalTag"
        done
    done
else
    echo "[Info] RUN_DATAVMC_PLOTS=$RUN_DATAVMC_PLOTS; skip 2_plot_dataVmc.py"
fi

if [[ "$RUN_OPTIMIZATION" == "1" ]]; then
    for finalTag in nominal sideband_rwgt; do
        python3 "$scriptsDir/ALP_Optimization.py" \
            -y run3 \
            -o "${outputDir}/optimize_run3UL_${finalTag}" \
            --region 1 \
            -p --sigVSscore -s --doOpt -c 1 \
            --inputTag "$finalTag"
    done
else
    echo "[Info] RUN_OPTIMIZATION=$RUN_OPTIMIZATION; skip ALP_Optimization.py"
fi
