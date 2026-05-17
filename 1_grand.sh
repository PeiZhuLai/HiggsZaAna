#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"

echo_step() {
    echo
    echo "==== $* ===="
}

has_queue_rows() {
    local joblist="$1"
    [[ -s "$joblist" ]] || return 1
    awk 'NF && $1 !~ /^#/ { found=1; exit } END { exit(found ? 0 : 1) }' "$joblist"
}

submit_and_wait() {
    local submit_file="$1"
    local log_pattern="$2"
    local submit_output cluster_id log_file

    mkdir -p "$(dirname "$(printf "$log_pattern" "0")")"
    submit_output="$(condor_submit -terse "$submit_file")"
    cluster_id="$(sed -E 's/^[^0-9]*([0-9]+).*/\1/' <<< "$submit_output")"
    log_file="$(printf "$log_pattern" "$cluster_id")"
    echo "[Condor] submitted cluster ${cluster_id}: ${submit_file}"
    condor_wait "$log_file"
    echo "[Condor] cluster ${cluster_id} finished"
}

# The story after HDNA parquets have been produced.

# Environmnet setup
# use-anaconda
# anaconda
# conda activate higgs-alp-ana

# p2root
echo_step "HiggsZaAna: p2root"
cd "${PROJECT_DIR}/Parquet2Rootfile"
bash 1_run_P2Root.sh
bash 2_prepare_rootfile.sh

# Train MVA
echo_step "HiggsZaAna: train MVA"
cd "${PROJECT_DIR}/HZaMVA/scripts"
python3 1_make_sideband_reweight.py
bash 2_train.sh
bash 3_save_model.sh

# p2root for MVA Score
echo_step "HiggsZaAna: p2root for MVA score"
cd "${PROJECT_DIR}/Parquet2Rootfile/Condor"
python3 1_make_joblist.py

if has_queue_rows joblist.tsv; then
    submit_and_wait "2_submit.sub" "${PROJECT_DIR}/Parquet2Rootfile/Condor/logs/%s.log"
else
    echo "[Condor] no p2root jobs to submit"
fi

python3 3_condor_resubmit.py --no-submit
if has_queue_rows joblist.tsv; then
    submit_and_wait "2_submit.sub" "${PROJECT_DIR}/Parquet2Rootfile/Condor/logs/%s.log"
else
    echo "[Condor] no p2root resubmit jobs to submit"
fi

bash 4_prepaare_2024DYJetsToLL.sh

# Determine MVA Cut
echo_step "HiggsZaAna: determine MVA cut"
cd "${PROJECT_DIR}/Plot/Condor"

NO_SUBMIT=1 bash 1_submit_dataVmc_condor.sh
if has_queue_rows dataVmc_jobs.txt; then
    submit_and_wait "dataVmc.submit" "${PROJECT_DIR}/Plot/Condor/logs/dataVmc.%s.log"
else
    echo "[Condor] no dataVmc jobs to submit"
fi

NO_SUBMIT=1 bash 2_resubmit_dataVmc_condor.sh
if has_queue_rows dataVmc_resubmit_jobs.txt; then
    submit_and_wait "dataVmc_resubmit.submit" "${PROJECT_DIR}/Plot/Condor/logs/dataVmc.%s.log"
else
    echo "[Condor] no dataVmc resubmit jobs to submit"
fi

bash 3_merge_dataVmc_condor.sh

# Plot associated plots
echo_step "HiggsZaAna: plot associated plots"
cd "${PROJECT_DIR}/Plot"
bash 1_runPlot.sh

echo_step "Switch to flashggFinalFit environment"
conda deactivate
cd "${PROJECT_DIR}"

####--------------------------------------------
####----------- flashggFinalFit ----------------
####--------------------------------------------

baseDir=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd "$baseDir"
cmsenv

### ----------- MVA Cut 
echo_step "flashggFinalFit: MVA cut"
python3 "$baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py" &
bash "$baseDir/shellScripts/mva/run_apply_bdt_sig_6jobs.sh"

# wait
# ### ----------- Tree2WS
echo_step "flashggFinalFit: Tree2WS"
cd "$baseDir/Trees2WS"
bash "$baseDir/Trees2WS/run_tree2ws.sh"

### ----------- Background 
# cd $baseDir/shellScripts
# sh $baseDir/shellScripts/bkg/fit_bkg.sh
### ----------- Background Condor
# cd $baseDir/shellScripts
# bash $baseDir/shellScripts/bkg/Condor/subjob_bkg.sh
echo_step "flashggFinalFit: background"
bash "$baseDir/shellScripts/bkg/Condor/collect_bkg_results.sh"
### ----------- Background Condor run locally
# cd $baseDir/shellScripts
# bash bkg/Condor/subjob_bkg.sh
# bash bkg/Condor/collect_bkg_results.sh

### ----------- Signal
echo_step "flashggFinalFit: signal"
cd "$baseDir/shellScripts"
bash "$baseDir/shellScripts/sig_sys/1_runjob_sig_fTest.sh"
bash "$baseDir/shellScripts/sig_sys/2_runjob_sig_calcPhotonSyst.sh"
bash "$baseDir/shellScripts/sig_sys/3_runjob_sig_signalFit.sh"
bash "$baseDir/shellScripts/sig_sys/4_runjob_sig_RunPlotter.sh"
bash "$baseDir/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh"

### ----------- Datacard
echo_step "flashggFinalFit: datacard"
cd "$baseDir/Datacard"
sh 1_runjob_gen_datacard_makeYields.sh
sh 2_runjob_gen_datacard_makeDatacard.sh
sh 3_rysn_datacard.sh

# ### ----------- Combine Limits
echo_step "flashggFinalFit: combine limits"
cd "$baseDir/Combine"
sh 1_makeLimits.sh
sh 2_text2ws.sh

# ### ----------- Plot Limts
echo_step "flashggFinalFit: plot limits"
cd "$baseDir/Plots"
sh 1_runLimitsPlot.sh

### ----------- Impact Plot
# cd $baseDir/Combine
# sh 3_expectedImpact.sh
### ----------- Impact Plot Condor (Generate ws for bias study)
echo_step "flashggFinalFit: impact plot condor"
bash "$baseDir/shellScripts/impact/Condor/subjob_expectedImpact.sh"

### ----------- Bias Study
# cd $baseDir/Combine/Checks/Bias_nominal
# sh 1_bias_study.sh
### ----------- Bias Study Condor (not working, need to check)
echo_step "flashggFinalFit: bias study condor"
bash "$baseDir/shellScripts/bias/Condor/subjob_bias_study.sh"

### ----------- Collect Bkg Fit Summary
echo_step "flashggFinalFit: collect background fit summary"
bash "$baseDir/shellScripts/bkg/Condor/collect_bkg_results.sh"

cd "$baseDir/shellScripts"


####--------------------------------------------
####-------------- Update AN -------------------
####--------------------------------------------

echo_step "Exit CMSSW environment"
if command -v scram >/dev/null 2>&1; then
    eval "$(scram unsetenv -sh)"
else
    echo "[Warning] scram not found; skip CMSSW environment cleanup"
fi

ANbashDir=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
echo_step "Update AN"
cd "$ANbashDir"
git pull

bash sync_figures.sh
bash compile.sh

git add .
git commit -m "AN and output"
git push
