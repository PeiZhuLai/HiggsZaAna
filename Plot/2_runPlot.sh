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
