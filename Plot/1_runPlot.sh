#!/bin/bash
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
scriptsDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts'
outputDir='/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots'
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/lib


# ## MVA Score signal region
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 1 --ln &

# ## MVA Score control region
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 2 --ln &

# wait

# ## ALP Optimization 2 categories
# python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 2

## ALP Optimization 1 category
python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 1

## Signal efficinecy after MVA cut
python3 $scriptsDir/collect_MVAcut_points_run3.py

## Batch 1
python3 $scriptsDir/signal_eff_sumw.py &
python $scriptsDir/BDT_ma_2D_lib.py &
python3 $scriptsDir/table_MVAScore_sigEff_bkgYield.py &

wait

## Batch 2
python3 $scriptsDir/plot_cutflowVmA.py &
python $scriptsDir/table_interpolate_bkgYield_1.py
if [ $? -eq 0 ]; then
    python $scriptsDir/table_interpolate_bkgYield_2.py &
else
    echo "[Batch 2] table_interpolate_bkgYield_1.py failed; skip table_interpolate_bkgYield_2.py"
fi
wait

## Batch 3
python3 $scriptsDir/plot_preselectSigEffVmA.py &
python3 $scriptsDir/plot_preselectSigEffSumwVmA.py &
python3 $scriptsDir/plot_MVASigEffVmA.py &

wait

## Batch 4
python3 $scriptsDir/plot_phidVmA.py &
python3 $scriptsDir/plot_eachphidVmA.py &
python3 $scriptsDir/plot_phidsigniVmA.py &

wait

## Batch 5
python3 $scriptsDir/plot_SIP3DsigniVmA.py &
python3 $scriptsDir/plot_dREff.py &
python3 $scriptsDir/plot_dREffBar.py &

wait

## Batch 6
python3 $scriptsDir/plot_trigEffVlepPt.py &
python3 $scriptsDir/plot_mAmigratedBar.py &
python3 $scriptsDir/plot_mAmigratedHist.py &
python3 $scriptsDir/plot_mAmigratedMatrix.py &
wait

python3 $scriptsDir/plot_bkgmcSculptingCheck.py &
python3 $scriptsDir/plot_mH_phopT_2D.py &
python3 $scriptsDir/plot_mH_mZ_2D.py &
wait

## Validate SF
python3 $scriptsDir/plot_SF_validation.py -y run3 -m --ln -b &

## fast variable plot
python3 $scriptsDir/plot_fast_variable_dataVmc.py -y run3 -m --ln -b --skip-sys &
wait

# ## MVA Score whole region
# python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --ln -b &
# wait
