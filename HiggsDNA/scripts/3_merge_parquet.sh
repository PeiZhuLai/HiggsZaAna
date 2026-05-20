#!/bin/bash

scriptsDir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts"

# Nominal Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Bkg_MC

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Sig_MC

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Data

# Control Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Bkg_MC_control

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Sig_MC_control

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Data_control
