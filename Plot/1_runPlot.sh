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

scriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts'
outputDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots'
variablesDir="${outputDir}/variables_dataVmc"
sidebandReweightJson='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json'
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib"

RUN_DATAVMC_PLOTS="${RUN_DATAVMC_PLOTS:-1}"
logDir="${outputDir}/logs_split"
mkdir -p "$logDir"

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

# 0: Full Region, 1: Signal Region, 2: Contral Region
# ## ALP Optimization 2 categories
python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 2 --inputTag sideband_rwgt

# ## ALP Optimization 1 category
python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 1 --inputTag sideband_rwgt

# ## Update cutflow
python3 /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/4_collect_cutflow.py
python3 /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts/5_merge_cutflow.py

# # Signal efficinecy after MVA cut
python3 $scriptsDir/collect_MVAcut_points_run3.py

## MVA model diagnostic plots (ROC, correlation matrices, per-feature distributions,
## feature importance, BDT score, and the low-mass overtraining plot). Regenerated from the
## DEPLOYED model pkls (NO retraining), so the AN figures always match the model that produced
## the limits. Scripts live in HZaMVA/scripts (must run from there for hza_features import).
mvaScriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts'
mvaUsingDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using'
mvaPlotDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/plots_MVA'
(
  cd "$mvaScriptsDir" || exit 1
  python3 plot_mva_diagnostics.py --model "$mvaUsingDir/model_Za_BDT_lowmass_run3.pkl"  --meta "$mvaUsingDir/model_Za_BDT_lowmass_run3.meta.json"  --outdir "$mvaPlotDir/run3_lowmass"
  python3 plot_mva_diagnostics.py --model "$mvaUsingDir/model_Za_BDT_highmass_run3.pkl" --meta "$mvaUsingDir/model_Za_BDT_highmass_run3.meta.json" --outdir "$mvaPlotDir/run3_highmass"
  python3 plot_overtrain_lowmass.py
)

## Batch 1
python3 $scriptsDir/signal_eff_sumw.py &
python $scriptsDir/BDT_ma_2D_lib.py &
python3 $scriptsDir/table_MVAScore_sigEff_bkgYield.py &
python3 $scriptsDir/plot_sigGenInfo.py &
wait

### Batch 2
python $scriptsDir/table_interpolate_bkgYield_1.py
python $scriptsDir/table_interpolate_bkgYield_2.py
wait

#### Batch 3
python3 $scriptsDir/plot_cutflowVmA.py &
python3 $scriptsDir/plot_preselectSigEffVmA.py &
python3 $scriptsDir/plot_preselectSigEffSumwVmA.py &
python3 $scriptsDir/plot_MVASigEffVmA.py &

wait

#### Batch 4
python3 $scriptsDir/plot_phidVmA.py &
python3 $scriptsDir/plot_eachphidVmA.py &
python3 $scriptsDir/plot_phidsigniVmA.py &

wait

#### Batch 5
python3 $scriptsDir/plot_SIP3DsigniVmA.py &
python3 $scriptsDir/plot_dREff.py &
python3 $scriptsDir/plot_dREffBar.py &

wait

### Batch 6
python3 $scriptsDir/plot_trigEffVlepPt.py &
python3 $scriptsDir/plot_mAmigratedBar.py &
python3 $scriptsDir/plot_mAmigratedHist.py &
python3 $scriptsDir/plot_mAmigratedMatrix.py &
wait

python3 $scriptsDir/plot_mH_phopT_2D.py &
python3 $scriptsDir/plot_mH_mZ_2D.py &
python3 $scriptsDir/plot_bkgmcSculptingCheck.py &
wait

### Validate SF
python3 $scriptsDir/plot_SF_validation.py -y run3 -m --ln -b &

### fast variable plot
python3 $scriptsDir/plot_fast_variable_dataVmc.py -y run3 -m --ln -b --skip-sys &
wait



# ##########################################################################################
# ## MVA Score whole region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --ln -b &
# wait

# ## MVA Score signal region
# # MVA Score signal region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --region 1 --ln &

# # MVA Score control region NFlow
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3_NFlow -m --region 2 --ln &

# ## ALP Optimization 2 categories NFlow
# python3 $scriptsDir/ALP_Optimization.py -y run3_NFlow -o $outputDir/optimize_run3UL_NFlow --region 2 -p --sigVSscore -s --doOpt -c 2

# ## ALP Optimization 1 category NFlow
# python3 $scriptsDir/ALP_Optimization.py -y run3_NFlow -o $outputDir/optimize_run3UL_NFlow --region 2 -p --sigVSscore -s --doOpt -c 1
