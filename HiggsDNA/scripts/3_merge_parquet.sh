#!/bin/bash

scriptsDir="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts"

# Nominal Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Bkg_MC --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Sig_MC --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp/Data --force

# Control Samples

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Bkg_MC_control --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Sig_MC_control --force

PYTHONPATH=$PWD python $scriptsDir/A_merge_parquet_outputs.py /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_control_tmp/Data_control --force
