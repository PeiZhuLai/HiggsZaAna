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


## MVA Score whole region
python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --ln -b

## MVA Score signal region
python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 1 --ln

## MVA Score control region
python3 $scriptsDir/plot_variable_dataVmc.py -y run3 -m --region 2 --ln

## ALP Optimization 2 categories
python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 2

## ALP Optimization 1 category
python3 $scriptsDir/ALP_Optimization.py -y run3 -o $outputDir/optimize_run3UL --region 2 -p --sigVSscore -s --doOpt -c 1

## Signal efficinecy after MVA cut
python3 $scriptsDir/collect_MVAcut_points_run3.py
python3 $scriptsDir/signal_eff_sumw.py

## MVA ScoreVmA
python $scriptsDir/BDT_ma_2D_lib.py

## MVAScore signal efficiency and background yield table
python3 $scriptsDir/table_MVAScore_sigEff_bkgYield.py

## Interpolate background yield table
python $scriptsDir/table_interpolate_bkgYield_1.py
python $scriptsDir/table_interpolate_bkgYield_2.py

## Plot cutflow vs ALP mass
python3 $scriptsDir/plot_cutflowVmA.py

## Unweight preselection efficiency
python3 $scriptsDir/plot_preselectSigEffVmA.py

## Weighted preselection efficiency
python3 $scriptsDir/plot_preselectSigEffSumwVmA.py

## Plot MVA cut efficiency vs ALP mass
python3 $scriptsDir/plot_MVASigEffVmA.py

## Photon ID signal efficiency and significance VS mA 
python3 $scriptsDir/plot_phidVmA.py
python3 $scriptsDir/plot_eachphidVmA.py
python3 $scriptsDir/plot_phidsigniVmA.py

## electron SIP3D study (eff_sig / eff_bkg) vs mA || Signficance vs mA 
python3 $scriptsDir/plot_SIP3DsigniVmA.py

## dR efficiency and dR efficiency bar plot
python3 $scriptsDir/plot_dREff.py
python3 $scriptsDir/plot_dREffBar.py

## Plot trigger efficiency vs mA
python3 $scriptsDir/plot_trigEffVmA.py

## Plot MVA score distribution for different mA
python3 $scriptsDir/plot_mAmigratedBar.py
#----------------------------------------------------------------------------
