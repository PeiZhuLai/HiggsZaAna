#!/usr/bin/env bash
# Data/MC plots of the selected BDT inputs (plot_datamc_inputs.py), env higgs-alp-ana activated (PyROOT needs the
# env compiler on PATH). Usage: run_datamc.sh <regions>   e.g. sideband | zmmg | sideband,zmmg
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
python plot_datamc_inputs.py --regions "$1"
