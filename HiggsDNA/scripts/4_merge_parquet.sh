#!/bin/bash

scriptsDir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts"
outputDir="/eos/project/h/htozg-dy-privatemc/pelai/HZa"

# Nominal Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA/Bkg_MC --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA/Sig_MC --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA/Data --force

# Control Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA_control/Bkg_MC_control --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA_control/Sig_MC_control --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py $outputDir/parquet_DNA_control/Data_control --force
