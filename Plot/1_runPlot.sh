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
#------------------------------Dry Run---------------------------------------
#1
# python scripts/plot_variable_dataVmc.py -y run3 --ln -b
#2
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --ln -b
#3
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --region 1 --ln
#4
# python3 scripts/plot_variable_dataVmc.py -y run3 -m --region 2 --ln
#5
# python3 scripts/ALP_Optimization.py -y run3 -o ./optimize_run3UL --region 0 -p --sigVSscore -s --doOpt -c 1
#----------------------------------------------------------------------------

# 0: Full Region, 1: Signal Region, 2: Contral Region
#-----------------------------Parallel Run-----------------------------------
set -euo pipefail

scriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts'
outputDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots'
variablesDir="${outputDir}/variables_dataVmc"
sidebandReweightJson='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json'
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib"

MAX_PLOT_JOBS="${MAX_PLOT_JOBS:-6}"
RUN_OPTIMIZATION="${RUN_OPTIMIZATION:-1}"
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
    local cmd=(python3 "$scriptsDir/plot_variable_dataVmc.py" -y run3 -m --ln --histOnly --samples "$samples" --outputTag "$partial_tag")

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

for finalTag in nominal sideband_rwgt; do
    for regionKey in SR CR mva; do
        merge_plot_output "$regionKey" "$finalTag"
    done
done

if [[ "$RUN_OPTIMIZATION" == "1" ]]; then
    for finalTag in nominal sideband_rwgt; do
        python3 "$scriptsDir/ALP_Optimization.py" \
            -y run3 \
            -o "${outputDir}/optimize_run3UL_${finalTag}" \
            --region 2 \
            -p --sigVSscore -s --doOpt -c 2 \
            --inputTag "$finalTag"
    done
else
    echo "[Info] RUN_OPTIMIZATION=$RUN_OPTIMIZATION; skip ALP_Optimization.py"
fi

#### Test Run with maxEvents
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 1 --ln --outputTag sideband_rwgt --useSidebandReweight --sidebandReweightJson $sidebandReweightJson  --maxEvents 100 &
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 2 --ln --outputTag sideband_rwgt --useSidebandReweight --sidebandReweightJson $sidebandReweightJson  --maxEvents 100 &
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --ln -b --outputTag sideband_rwgt --useSidebandReweight --sidebandReweightJson $sidebandReweightJson  --maxEvents 100 &
# wait

# ## ALP Optimization 2 categories
# RUN_OPTIMIZATION=1 bash /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/1_runPlot.sh

# ## ALP Optimization 1 category
# python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 1 --inputTag sideband_rwgt

# # Signal efficinecy after MVA cut
# python3 $scriptsDir/collect_MVAcut_points_run3.py

# ## Batch 1
# python3 $scriptsDir/signal_eff_sumw.py &
# python $scriptsDir/BDT_ma_2D_lib.py &
# python3 $scriptsDir/table_MVAScore_sigEff_bkgYield.py &
# python3 $scriptsDir/plot_sigGenInfo.py &
# wait

# ### Batch 2
# python $scriptsDir/table_interpolate_bkgYield_1.py
# python $scriptsDir/table_interpolate_bkgYield_2.py
# wait

# #### Batch 3
# python3 $scriptsDir/plot_cutflowVmA.py &
# python3 $scriptsDir/plot_preselectSigEffVmA.py &
# python3 $scriptsDir/plot_preselectSigEffSumwVmA.py &
# python3 $scriptsDir/plot_MVASigEffVmA.py &

# wait

# #### Batch 4
# python3 $scriptsDir/plot_phidVmA.py &
# python3 $scriptsDir/plot_eachphidVmA.py &
# python3 $scriptsDir/plot_phidsigniVmA.py &

# wait

# #### Batch 5
# python3 $scriptsDir/plot_SIP3DsigniVmA.py &
# python3 $scriptsDir/plot_dREff.py &
# python3 $scriptsDir/plot_dREffBar.py &

# wait

# ### Batch 6
# python3 $scriptsDir/plot_trigEffVlepPt.py &
# python3 $scriptsDir/plot_mAmigratedBar.py &
# python3 $scriptsDir/plot_mAmigratedHist.py &
# python3 $scriptsDir/plot_mAmigratedMatrix.py &
# wait

# python3 $scriptsDir/plot_mH_phopT_2D.py &
# python3 $scriptsDir/plot_mH_mZ_2D.py &
# python3 $scriptsDir/plot_bkgmcSculptingCheck.py &
# wait

# ### Validate SF
# python3 $scriptsDir/plot_SF_validation.py -y run3 -m --ln -b &

# ### fast variable plot
# python3 $scriptsDir/plot_fast_variable_dataVmc.py -y run3 -m --ln -b --skip-sys &
# wait


###########################################################################################
# ## MVA Score whole region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --ln -b &
# wait

# ## MVA Score signal region
## MVA Score signal region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --region 1 --ln &

## MVA Score control region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --region 2 --ln &

# ## ALP Optimization 2 categories NFlow
# python3 $scriptsDir/ALP_Optimization.py -y run3_NFlow -o $outputDir/optimize_run3UL_NFlow --region 2 -p --sigVSscore -s --doOpt -c 2

# ## ALP Optimization 1 category NFlow
# python3 $scriptsDir/ALP_Optimization.py -y run3_NFlow -o $outputDir/optimize_run3UL_NFlow --region 2 -p --sigVSscore -s --doOpt -c 1
