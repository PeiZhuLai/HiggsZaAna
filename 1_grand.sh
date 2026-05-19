#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"

RUN_HIGGSZA_P2ROOT="${RUN_HIGGSZA_P2ROOT:-0}"
RUN_HIGGSZA_TRAIN_MVA="${RUN_HIGGSZA_TRAIN_MVA:-0}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE="${RUN_HIGGSZA_P2ROOT_MVA_SCORE:-0}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE_INITIAL="${RUN_HIGGSZA_P2ROOT_MVA_SCORE_INITIAL:-0}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE_RESUBMIT="${RUN_HIGGSZA_P2ROOT_MVA_SCORE_RESUBMIT:-0}"
RUN_HIGGSZA_DETERMINE_MVA_CUT="${RUN_HIGGSZA_DETERMINE_MVA_CUT:-0}"
RUN_HIGGSZA_DETERMINE_MVA_CUT_INITIAL="${RUN_HIGGSZA_DETERMINE_MVA_CUT_INITIAL:-0}"
RUN_HIGGSZA_DETERMINE_MVA_CUT_RESUBMIT="${RUN_HIGGSZA_DETERMINE_MVA_CUT_RESUBMIT:-0}"
RUN_HIGGSZA_PLOT="${RUN_HIGGSZA_PLOT:-0}"
RUN_FLASHGG_ENV="${RUN_FLASHGG_ENV:-0}"
RUN_FLASHGG_MVA_CUT="${RUN_FLASHGG_MVA_CUT:-0}"
RUN_FLASHGG_TREE2WS="${RUN_FLASHGG_TREE2WS:-0}"
RUN_FLASHGG_BACKGROUND="${RUN_FLASHGG_BACKGROUND:-0}"
RUN_FLASHGG_SIGNAL="${RUN_FLASHGG_SIGNAL:-0}"
RUN_FLASHGG_DATACARD="${RUN_FLASHGG_DATACARD:-1}"
RUN_FLASHGG_COMBINE_LIMITS="${RUN_FLASHGG_COMBINE_LIMITS:-1}"
RUN_FLASHGG_PLOT_LIMITS="${RUN_FLASHGG_PLOT_LIMITS:-1}"
RUN_FLASHGG_IMPACT="${RUN_FLASHGG_IMPACT:-1}"
RUN_FLASHGG_BIAS="${RUN_FLASHGG_BIAS:-1}"
RUN_FLASHGG_COLLECT_BKG="${RUN_FLASHGG_COLLECT_BKG:-1}"
RUN_EXIT_CMSSW_ENV="${RUN_EXIT_CMSSW_ENV:-1}"
RUN_UPDATE_AN="${RUN_UPDATE_AN:-1}"
P2ROOT_RESUBMIT_MAX_ATTEMPTS="${P2ROOT_RESUBMIT_MAX_ATTEMPTS:-5}"
P2ROOT_RESUBMIT_CHECK_ROOT="${P2ROOT_RESUBMIT_CHECK_ROOT:-0}"
DATA_VMC_RESUBMIT_MAX_ATTEMPTS="${DATA_VMC_RESUBMIT_MAX_ATTEMPTS:-5}"

rm -fr /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal
rm -fr /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal
rm -fr /eos/home-p/pelai/HZa/root_MVAcut
rm -fr /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor/logs
rm -fr /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/log_file
rm -fr /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc/*_part_*.root



echo_step() {
    echo
    echo "==== $* ===="
}

echo_skip() {
    echo_step "Skip: $*"
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

p2root_resubmit_until_done() {
    local attempt
    local resubmit_args=(--no-submit)

    if [[ "$P2ROOT_RESUBMIT_CHECK_ROOT" == "1" ]]; then
        resubmit_args+=(--check-root)
    fi

    for ((attempt = 1; attempt <= P2ROOT_RESUBMIT_MAX_ATTEMPTS; attempt++)); do
        echo_step "HiggsZaAna: p2root resubmit check ${attempt}/${P2ROOT_RESUBMIT_MAX_ATTEMPTS}"
        python3 3_condor_resubmit.py "${resubmit_args[@]}"

        if ! has_queue_rows joblist.tsv; then
            echo "[Condor] no p2root resubmit jobs to submit"
            return 0
        fi

        submit_and_wait "2_submit.sub" "${PROJECT_DIR}/Parquet2Rootfile/Condor/logs/%s.log"
    done

    echo_step "HiggsZaAna: final p2root resubmit check"
    python3 3_condor_resubmit.py "${resubmit_args[@]}"
    if has_queue_rows joblist.tsv; then
        echo "[ERROR] p2root still has unfinished outputs after ${P2ROOT_RESUBMIT_MAX_ATTEMPTS} resubmit attempts." >&2
        echo "[ERROR] Remaining jobs are listed in ${PROJECT_DIR}/Parquet2Rootfile/Condor/joblist.tsv" >&2
        return 1
    fi

    echo "[Condor] no p2root resubmit jobs to submit"
}

data_vmc_resubmit_until_done() {
    local attempt

    for ((attempt = 1; attempt <= DATA_VMC_RESUBMIT_MAX_ATTEMPTS; attempt++)); do
        echo_step "HiggsZaAna: dataVmc resubmit check ${attempt}/${DATA_VMC_RESUBMIT_MAX_ATTEMPTS}"
        NO_SUBMIT=1 bash 2_resubmit_dataVmc_condor.sh

        if ! has_queue_rows dataVmc_resubmit_jobs.txt; then
            echo "[Condor] no dataVmc resubmit jobs to submit"
            return 0
        fi

        submit_and_wait "dataVmc_resubmit.submit" "${PROJECT_DIR}/Plot/Condor/logs/dataVmc.%s.log"
    done

    echo_step "HiggsZaAna: final dataVmc resubmit check"
    NO_SUBMIT=1 bash 2_resubmit_dataVmc_condor.sh
    if has_queue_rows dataVmc_resubmit_jobs.txt; then
        echo "[ERROR] dataVmc still has unfinished outputs after ${DATA_VMC_RESUBMIT_MAX_ATTEMPTS} resubmit attempts." >&2
        echo "[ERROR] Remaining jobs are listed in ${PROJECT_DIR}/Plot/Condor/dataVmc_resubmit_jobs.txt" >&2
        echo "[ERROR] Missing outputs are listed in ${PROJECT_DIR}/Plot/Condor/dataVmc_missing_outputs.txt" >&2
        return 1
    fi

    echo "[Condor] no dataVmc resubmit jobs to submit"
}

# The story after HDNA parquets have been produced.

# Environmnet setup
# use-anaconda
# anaconda
# conda activate higgs-alp-ana

# p2root
if [[ "$RUN_HIGGSZA_P2ROOT" == "1" ]]; then
    echo_step "HiggsZaAna: p2root"
    cd "${PROJECT_DIR}/Parquet2Rootfile"
    bash 1_run_P2Root.sh
    bash 2_prepare_rootfile.sh
else
    echo_skip "HiggsZaAna: p2root"
fi

# Train MVA
if [[ "$RUN_HIGGSZA_TRAIN_MVA" == "1" ]]; then
    echo_step "HiggsZaAna: train MVA"
    cd "${PROJECT_DIR}/HZaMVA/scripts"
    python3 1_make_sideband_reweight.py
    bash 2_train.sh
    bash 3_save_model.sh
else
    echo_skip "HiggsZaAna: train MVA"
fi

# p2root for MVA Score
if [[ "$RUN_HIGGSZA_P2ROOT_MVA_SCORE" == "1" ]]; then
    echo_step "HiggsZaAna: p2root for MVA score"
    cd "${PROJECT_DIR}/Parquet2Rootfile/Condor"

    if [[ "$RUN_HIGGSZA_P2ROOT_MVA_SCORE_INITIAL" == "1" ]]; then
        echo_step "HiggsZaAna: p2root for MVA score initial submit"
        python3 1_make_joblist.py

        if has_queue_rows joblist.tsv; then
            submit_and_wait "2_submit.sub" "${PROJECT_DIR}/Parquet2Rootfile/Condor/logs/%s.log"
        else
            echo "[Condor] no p2root jobs to submit"
        fi
    else
        echo_skip "HiggsZaAna: p2root for MVA score initial submit"
    fi

    if [[ "$RUN_HIGGSZA_P2ROOT_MVA_SCORE_RESUBMIT" == "1" ]]; then
        p2root_resubmit_until_done
    else
        echo_skip "HiggsZaAna: p2root for MVA score resubmit"
    fi

    bash 4_prepaare_2024DYJetsToLL.sh
else
    echo_skip "HiggsZaAna: p2root for MVA score"
fi

# Determine MVA Cut
if [[ "$RUN_HIGGSZA_DETERMINE_MVA_CUT" == "1" ]]; then
    echo_step "HiggsZaAna: determine MVA cut"
    cd "${PROJECT_DIR}/Plot/Condor"

    if [[ "$RUN_HIGGSZA_DETERMINE_MVA_CUT_INITIAL" == "1" ]]; then
        echo_step "HiggsZaAna: determine MVA cut initial submit"
        NO_SUBMIT=1 bash 1_submit_dataVmc_condor.sh
        if has_queue_rows dataVmc_jobs.txt; then
            submit_and_wait "dataVmc.submit" "${PROJECT_DIR}/Plot/Condor/logs/dataVmc.%s.log"
        else
            echo "[Condor] no dataVmc jobs to submit"
        fi
    else
        echo_skip "HiggsZaAna: determine MVA cut initial submit"
    fi

    if [[ "$RUN_HIGGSZA_DETERMINE_MVA_CUT_RESUBMIT" == "1" ]]; then
        data_vmc_resubmit_until_done
    else
        echo_skip "HiggsZaAna: determine MVA cut resubmit"
    fi

    bash 3_merge_dataVmc_condor.sh
else
    echo_skip "HiggsZaAna: determine MVA cut"
fi

# Plot associated plots
if [[ "$RUN_HIGGSZA_PLOT" == "1" ]]; then
    echo_step "HiggsZaAna: plot associated plots"
    cd "${PROJECT_DIR}/Plot"
    bash 1_runPlot.sh
else
    echo_skip "HiggsZaAna: plot associated plots"
fi

if [[ "$RUN_FLASHGG_ENV" == "1" ]]; then
    echo_step "Switch to flashggFinalFit environment"
    conda deactivate
    cd "${PROJECT_DIR}"
else
    echo_skip "Switch to flashggFinalFit environment"
fi

####--------------------------------------------
####----------- flashggFinalFit ----------------
####--------------------------------------------

baseDir=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
if [[ "$RUN_FLASHGG_ENV" == "1" ]]; then
    cd "$baseDir"
    cmsenv
fi

### ----------- MVA Cut 
if [[ "$RUN_FLASHGG_MVA_CUT" == "1" ]]; then
    echo_step "flashggFinalFit: MVA cut"
    python3 "$baseDir/MVAcut/run3_ReReco_Sys/scripts/apply_bdt_data.py" &
    bash "$baseDir/shellScripts/mva/run_apply_bdt_sig_6jobs.sh"
else
    echo_skip "flashggFinalFit: MVA cut"
fi

# wait
# ### ----------- Tree2WS
if [[ "$RUN_FLASHGG_TREE2WS" == "1" ]]; then
    echo_step "flashggFinalFit: Tree2WS"
    cd "$baseDir/Trees2WS"
    bash "$baseDir/Trees2WS/run_tree2ws.sh"
else
    echo_skip "flashggFinalFit: Tree2WS"
fi

### ----------- Background 
# cd $baseDir/shellScripts
# sh $baseDir/shellScripts/bkg/fit_bkg.sh
### ----------- Background Condor
# cd $baseDir/shellScripts
# bash $baseDir/shellScripts/bkg/Condor/subjob_bkg.sh
if [[ "$RUN_FLASHGG_BACKGROUND" == "1" ]]; then
    echo_step "flashggFinalFit: background"
    bash "$baseDir/shellScripts/bkg/Condor/subjob_bkg.sh"
    bash "$baseDir/shellScripts/bkg/Condor/collect_bkg_results.sh"
else
    echo_skip "flashggFinalFit: background"
fi
### ----------- Background Condor run locally
# cd $baseDir/shellScripts
# bash bkg/Condor/subjob_bkg.sh
# bash bkg/Condor/collect_bkg_results.sh

### ----------- Signal
if [[ "$RUN_FLASHGG_SIGNAL" == "1" ]]; then
    echo_step "flashggFinalFit: signal"
    cd "$baseDir/shellScripts"
    bash "$baseDir/shellScripts/sig_sys/1_runjob_sig_fTest.sh"
    bash "$baseDir/shellScripts/sig_sys/2_runjob_sig_calcPhotonSyst.sh"
    bash "$baseDir/shellScripts/sig_sys/3_runjob_sig_signalFit.sh"
    bash "$baseDir/shellScripts/sig_sys/4_runjob_sig_RunPlotter.sh"
    bash "$baseDir/shellScripts/sig_sys/5_runjob_sig_plotEffSigma.sh"
else
    echo_skip "flashggFinalFit: signal"
fi

### ----------- Datacard
if [[ "$RUN_FLASHGG_DATACARD" == "1" ]]; then
    echo_step "flashggFinalFit: datacard"
    cd "$baseDir/Datacard"
    sh 1_runjob_gen_datacard_makeYields.sh
    sh 2_runjob_gen_datacard_makeDatacard.sh
    sh 3_rysn_datacard.sh
else
    echo_skip "flashggFinalFit: datacard"
fi

# ### ----------- Combine Limits
if [[ "$RUN_FLASHGG_COMBINE_LIMITS" == "1" ]]; then
    echo_step "flashggFinalFit: combine limits"
    cd "$baseDir/Combine"
    sh 1_makeLimits.sh
    sh 2_text2ws.sh
else
    echo_skip "flashggFinalFit: combine limits"
fi

# ### ----------- Plot Limts
if [[ "$RUN_FLASHGG_PLOT_LIMITS" == "1" ]]; then
    echo_step "flashggFinalFit: plot limits"
    cd "$baseDir/Plots"
    sh 1_runLimitsPlot.sh
else
    echo_skip "flashggFinalFit: plot limits"
fi

### ----------- Impact Plot
# cd $baseDir/Combine
# sh 3_expectedImpact.sh
### ----------- Impact Plot Condor (Generate ws for bias study)
if [[ "$RUN_FLASHGG_IMPACT" == "1" ]]; then
    echo_step "flashggFinalFit: impact plot condor"
    bash "$baseDir/shellScripts/impact/Condor/subjob_expectedImpact.sh"
else
    echo_skip "flashggFinalFit: impact plot condor"
fi

### ----------- Bias Study
# cd $baseDir/Combine/Checks/Bias_nominal
# sh 1_bias_study.sh
### ----------- Bias Study Condor (not working, need to check)
if [[ "$RUN_FLASHGG_BIAS" == "1" ]]; then
    echo_step "flashggFinalFit: bias study condor"
    bash "$baseDir/shellScripts/bias/Condor/subjob_bias_study.sh"
else
    echo_skip "flashggFinalFit: bias study condor"
fi

### ----------- Collect Bkg Fit Summary
if [[ "$RUN_FLASHGG_COLLECT_BKG" == "1" ]]; then
    echo_step "flashggFinalFit: collect background fit summary"
    bash "$baseDir/shellScripts/bkg/Condor/collect_bkg_results.sh"
else
    echo_skip "flashggFinalFit: collect background fit summary"
fi

if [[ "$RUN_FLASHGG_ENV" == "1" ]]; then
    cd "$baseDir/shellScripts"
fi


####--------------------------------------------
####-------------- Update AN -------------------
####--------------------------------------------

if [[ "$RUN_EXIT_CMSSW_ENV" == "1" ]]; then
    echo_step "Exit CMSSW environment"
    if command -v scram >/dev/null 2>&1; then
        eval "$(scram unsetenv -sh)"
    else
        echo "[Warning] scram not found; skip CMSSW environment cleanup"
    fi
else
    echo_skip "Exit CMSSW environment"
fi

ANbaseDir=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
if [[ "$RUN_UPDATE_AN" == "1" ]]; then
    echo_step "Update AN"
    cd "$ANbaseDir"
    git pull

    bash sync_figures.sh
    bash compile.sh

    git add .
    git commit -m "AN and output"
    git push
else
    echo_skip "Update AN"
fi

cd "$PROJECT_DIR"
