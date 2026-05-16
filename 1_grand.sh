#!/bin/bash                                                                                                                                                                       
# The story after HDNA parquets have been produced.

# Environmnet setup
# use-anaconda
# anaconda
# conda activate higgs-alp-ana

# p2root
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
bash 1_run_P2Root.sh
bash 2_prepare_rootfile.sh

# Train MVA
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
python3 1_make_sideband_reweight.py
bash 2_train.sh
bash 3_save_model.sh

# p2root for MVA Score
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
python3 1_make_joblist.py
condor_submit 2_submit.sub
python3 3_condor_resubmit.py 
bash 4_prepaare_2024DYJetsToLL.sh

# Determine MVA Cut
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/Condor

bash 1_submit_dataVmc_condor.sh
bash 2_resubmit_dataVmc_condor.sh
bash 3_merge_dataVmc_condor.sh

# Plot associated plots
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
bash 1_runPlot.sh