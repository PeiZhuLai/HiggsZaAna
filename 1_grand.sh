#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
HIGGSZA_EOS_DIR="${HIGGSZA_EOS_DIR:-/eos/home-p/pelai/HZa}"
ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
ANACONDA_CONDA_SH="${ANACONDA_CONDA_SH:-/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh}"
MINICONDA_SETUP="${MINICONDA_SETUP:-}"
MINICONDA_CONDA_SH="${MINICONDA_CONDA_SH:-}"

RUN_HIGGSZA_CLEAN_OUTPUTS="${RUN_HIGGSZA_CLEAN_OUTPUTS:-0}"
RUN_HIGGSZA_MERGE_PARQUET="${RUN_HIGGSZA_MERGE_PARQUET:-0}"
RUN_HIGGSZA_P2ROOT="${RUN_HIGGSZA_P2ROOT:-0}"
RUN_HIGGSZA_TRAIN_MVA="${RUN_HIGGSZA_TRAIN_MVA:-1}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE="${RUN_HIGGSZA_P2ROOT_MVA_SCORE:-1}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE_INITIAL="${RUN_HIGGSZA_P2ROOT_MVA_SCORE_INITIAL:-1}"
RUN_HIGGSZA_P2ROOT_MVA_SCORE_RESUBMIT="${RUN_HIGGSZA_P2ROOT_MVA_SCORE_RESUBMIT:-1}"
RUN_HIGGSZA_DETERMINE_MVA_CUT="${RUN_HIGGSZA_DETERMINE_MVA_CUT:-1}"
RUN_HIGGSZA_DETERMINE_MVA_CUT_INITIAL="${RUN_HIGGSZA_DETERMINE_MVA_CUT_INITIAL:-1}"
RUN_HIGGSZA_DETERMINE_MVA_CUT_RESUBMIT="${RUN_HIGGSZA_DETERMINE_MVA_CUT_RESUBMIT:-1}"
RUN_HIGGSZA_PLOT="${RUN_HIGGSZA_PLOT:-1}"
RUN_FLASHGG_ENV="${RUN_FLASHGG_ENV:-1}"
RUN_FLASHGG_MVA_CUT="${RUN_FLASHGG_MVA_CUT:-1}"
RUN_FLASHGG_TREE2WS="${RUN_FLASHGG_TREE2WS:-1}"
RUN_FLASHGG_BACKGROUND="${RUN_FLASHGG_BACKGROUND:-1}"
RUN_FLASHGG_SIGNAL="${RUN_FLASHGG_SIGNAL:-1}"
RUN_FLASHGG_DATACARD="${RUN_FLASHGG_DATACARD:-1}"
RUN_FLASHGG_COMBINE_LIMITS="${RUN_FLASHGG_COMBINE_LIMITS:-1}"
RUN_FLASHGG_PLOT_LIMITS="${RUN_FLASHGG_PLOT_LIMITS:-1}"
RUN_FLASHGG_IMPACT="${RUN_FLASHGG_IMPACT:-1}"
RUN_FLASHGG_BIAS="${RUN_FLASHGG_BIAS:-1}"
RUN_FLASHGG_COLLECT_BKG="${RUN_FLASHGG_COLLECT_BKG:-1}"
RUN_EXIT_CMSSW_ENV="${RUN_EXIT_CMSSW_ENV:-1}"
RUN_UPDATE_AN="${RUN_UPDATE_AN:-1}"
P2ROOT_RESUBMIT_MAX_ATTEMPTS="${P2ROOT_RESUBMIT_MAX_ATTEMPTS:-5}"
P2ROOT_RESUBMIT_CHECK_ROOT="${P2ROOT_RESUBMIT_CHECK_ROOT:-5}"
DATA_VMC_RESUBMIT_MAX_ATTEMPTS="${DATA_VMC_RESUBMIT_MAX_ATTEMPTS:-5}"

echo_step() {
    echo
    echo "==== $* ===="
}

echo_skip() {
    echo_step "Skip: $*"
}

init_conda_shell() {
    local candidate conda_hook
    local candidates=("$@")

    if declare -F conda >/dev/null 2>&1; then
        return 0
    fi

    if command -v conda >/dev/null 2>&1; then
        set +u
        if conda_hook="$(conda shell.bash hook 2>/dev/null)"; then
            eval "$conda_hook"
        fi
        set -u
        if declare -F conda >/dev/null 2>&1; then
            return 0
        fi
    fi

    for candidate in "${candidates[@]}"; do
        [[ -n "$candidate" ]] || continue
        [[ -r "$candidate" ]] || continue

        set +u
        # shellcheck disable=SC1090
        if ! source "$candidate"; then
            set -u
            continue
        fi
        set -u

        if declare -F conda >/dev/null 2>&1; then
            return 0
        fi

        if command -v conda >/dev/null 2>&1; then
            set +u
            if conda_hook="$(conda shell.bash hook 2>/dev/null)"; then
                eval "$conda_hook"
            fi
            set -u
            if declare -F conda >/dev/null 2>&1; then
                return 0
            fi
        fi
    done

    return 1
}

activate_conda_env() {
    local env_name="$1"
    shift

    if [[ "${CONDA_DEFAULT_ENV:-}" == "$env_name" ]]; then
        echo "[ENV] conda env already active: ${env_name}"
        return 0
    fi

    if ! init_conda_shell "$@"; then
        echo "[ERROR] Failed to initialize conda shell for env '${env_name}'." >&2
        echo "[ERROR] Set ANACONDA_SETUP/MINICONDA_SETUP or *_CONDA_SH to a readable setup script." >&2
        return 1
    fi

    set +u
    conda activate "$env_name"
    set -u

    if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
        echo "[ERROR] Failed to activate conda env '${env_name}'." >&2
        return 1
    fi

    echo "[ENV] active CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV}"
    echo "[ENV] active CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
}

deactivate_conda_env_if_active() {
    if ! init_conda_shell \
        "$ANACONDA_SETUP" \
        "$ANACONDA_CONDA_SH" \
        "$MINICONDA_SETUP" \
        "$MINICONDA_CONDA_SH" \
        "${CONDA_EXE:+${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh}" \
        "${HOME}/anaconda3/etc/profile.d/conda.sh" \
        "${HOME}/miniconda3/etc/profile.d/conda.sh" \
        "/usr/etc/profile.d/conda.sh"; then
        return 0
    fi

    if [[ -z "${CONDA_DEFAULT_ENV:-}" && -z "${CONDA_PREFIX:-}" ]]; then
        return 0
    fi

    set +u
    conda deactivate || true
    set -u
}

has_queue_rows() {
    local joblist="$1"
    [[ -s "$joblist" ]] || return 1
    awk 'NF && $1 !~ /^#/ { found=1; exit } END { exit(found ? 0 : 1) }' "$joblist"
}

condor_log_cluster_status() {
    local cluster_id="$1"
    local log_file="$2"

    [[ -s "$log_file" ]] || return 1

    awk -v cluster_id="$cluster_id" '
        match($0, /^(000|001|005|009|012) \(([0-9]+)\.([0-9]+)\.000\)/, m) {
            if (m[2] != cluster_id) {
                next
            }

            code = m[1]
            proc = m[3]
            seen[proc] = 1

            if (code == "000") {
                status[proc] = "idle"
            } else if (code == "001") {
                status[proc] = "running"
            } else if (code == "005") {
                status[proc] = "completed"
            } else if (code == "009") {
                status[proc] = "aborted"
            } else if (code == "012") {
                status[proc] = "held"
            }
        }
        END {
            total = idle = running = completed = aborted = held = unknown = 0

            for (proc in seen) {
                total += 1
                if (status[proc] == "idle") {
                    idle += 1
                } else if (status[proc] == "running") {
                    running += 1
                } else if (status[proc] == "completed") {
                    completed += 1
                } else if (status[proc] == "aborted") {
                    aborted += 1
                } else if (status[proc] == "held") {
                    held += 1
                } else {
                    unknown += 1
                }
            }

            if (total == 0) {
                exit 1
            }

            printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", total, idle, running, completed, aborted, held, unknown
        }
    ' "$log_file"
}

wait_for_cluster_via_log() {
    local cluster_id="$1"
    local log_file="$2"
    local poll_seconds="${CONDOR_Q_POLL_SECONDS:-30}"
    local max_missing_before_seen="${CONDOR_LOG_MAX_MISSING_BEFORE_SEEN:-24}"
    local seen_cluster=0
    local missing_before_seen=0
    local status_line total idle running completed aborted held unknown

    echo "[Condor] fallback to user log polling for cluster ${cluster_id}"
    echo "[Condor] user log path: ${log_file}"

    while true; do
        if ! status_line="$(condor_log_cluster_status "$cluster_id" "$log_file" 2>/dev/null)"; then
            if [[ "$seen_cluster" -eq 1 ]]; then
                echo "[Condor][log] cluster ${cluster_id} has no newer parseable events; assuming it already reached a terminal state"
                return 0
            fi

            missing_before_seen=$((missing_before_seen + 1))
            echo "[Condor][log] cluster ${cluster_id} not visible in user log yet (${missing_before_seen}/${max_missing_before_seen})"

            if [[ "$missing_before_seen" -ge "$max_missing_before_seen" ]]; then
                echo "[WARN] cluster ${cluster_id} never became visible in user log; assuming it already left the queue." >&2
                return 0
            fi

            sleep 5
            continue
        fi

        IFS=$'\t' read -r total idle running completed aborted held unknown <<< "$status_line"
        seen_cluster=1
        missing_before_seen=0

        echo "[Condor][log] cluster ${cluster_id}: total=${total} idle=${idle} running=${running} completed=${completed} aborted=${aborted} held=${held} unknown=${unknown}"

        if [[ "$held" -gt 0 ]]; then
            echo "[ERROR] cluster ${cluster_id} has held job(s) according to user log ${log_file}." >&2
            return 1
        fi

        if [[ "$idle" -eq 0 && "$running" -eq 0 ]]; then
            echo "[Condor][log] cluster ${cluster_id} reached terminal state via user log"
            return 0
        fi

        sleep "$poll_seconds"
    done
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
    wait_for_cluster_via_queue "$cluster_id" "$log_file"
    echo "[Condor] cluster ${cluster_id} finished"
}

wait_for_cluster_via_queue() {
    local cluster_id="$1"
    local log_file="$2"
    # Poll the queue directly because condor_wait can fail on user-log parsing.
    local poll_seconds="${CONDOR_Q_POLL_SECONDS:-30}"
    local max_query_failures="${CONDOR_Q_MAX_FAILURES:-5}"
    local max_missing_before_seen="${CONDOR_Q_MAX_MISSING_BEFORE_SEEN:-12}"
    local query_failures=0
    local missing_before_seen=0
    local seen_cluster=0
    local q_output q_status
    local proc_id job_status hold_reason
    local total idle running held removing completed suspended unknown
    local first_hold_reason status_line

    echo "[Condor] waiting for cluster ${cluster_id} via condor_q polling"
    echo "[Condor] user log path: ${log_file}"

    while true; do
        set +e
        q_output="$(
            condor_q "$cluster_id" \
                -format '%d\t' ProcId \
                -format '%d\t' JobStatus \
                -format '%s\n' 'ifThenElse(isUndefined(HoldReason),"",HoldReason)'
        )"
        q_status=$?
        set -e

        if [[ "$q_status" -ne 0 ]]; then
            query_failures=$((query_failures + 1))
            echo "[WARN] condor_q failed for cluster ${cluster_id} (${query_failures}/${max_query_failures}); retry in ${poll_seconds}s" >&2
            if [[ -n "${q_output:-}" ]]; then
                echo "$q_output" >&2
            fi

            if [[ "$query_failures" -ge "$max_query_failures" ]]; then
                echo "[WARN] condor_q kept failing while waiting for cluster ${cluster_id}; falling back to user log polling." >&2
                wait_for_cluster_via_log "$cluster_id" "$log_file"
                return $?
            fi

            sleep "$poll_seconds"
            continue
        fi

        query_failures=0

        if [[ -z "${q_output//[[:space:]]/}" ]]; then
            if [[ "$seen_cluster" -eq 1 ]]; then
                echo "[Condor] cluster ${cluster_id} no longer in condor_q"
                return 0
            fi

            missing_before_seen=$((missing_before_seen + 1))
            echo "[Condor] cluster ${cluster_id} not visible in condor_q yet (${missing_before_seen}/${max_missing_before_seen})"

            if [[ "$missing_before_seen" -ge "$max_missing_before_seen" ]]; then
                echo "[WARN] cluster ${cluster_id} never became visible in condor_q; assuming it already left the queue." >&2
                return 0
            fi

            sleep 5
            continue
        fi

        seen_cluster=1
        missing_before_seen=0
        total=0
        idle=0
        running=0
        held=0
        removing=0
        completed=0
        suspended=0
        unknown=0
        first_hold_reason=""

        while IFS=$'\t' read -r proc_id job_status hold_reason; do
            [[ -n "${proc_id:-}" ]] || continue
            total=$((total + 1))

            case "$job_status" in
                1) idle=$((idle + 1)) ;;
                2|6) running=$((running + 1)) ;;
                3) removing=$((removing + 1)) ;;
                4) completed=$((completed + 1)) ;;
                5)
                    held=$((held + 1))
                    if [[ -z "$first_hold_reason" ]]; then
                        first_hold_reason="proc ${proc_id}: ${hold_reason:-<no HoldReason>}"
                    fi
                    ;;
                7) suspended=$((suspended + 1)) ;;
                *) unknown=$((unknown + 1)) ;;
            esac
        done <<< "$q_output"

        status_line="[Condor] cluster ${cluster_id}: total=${total} idle=${idle} running=${running} held=${held} removing=${removing} completed=${completed} suspended=${suspended} unknown=${unknown}"
        echo "$status_line"

        if [[ "$held" -gt 0 ]]; then
            echo "[ERROR] cluster ${cluster_id} has held job(s): ${first_hold_reason}" >&2
            return 1
        fi

        sleep "$poll_seconds"
    done
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

if [[ "$RUN_HIGGSZA_CLEAN_OUTPUTS" == "1" ]]; then
    echo_step "HiggsZaAna: clean previous outputs"
    rm -fr "${HIGGSZA_EOS_DIR}/root_P2Root/run3_bdt_inputs_nominal"
    rm -fr "${HIGGSZA_EOS_DIR}/root_P2Root/run3_bdt_scored_nominal"
    rm -fr "${HIGGSZA_EOS_DIR}/root_MVAcut"
    rm -fr "${PROJECT_DIR}/Parquet2Rootfile/Condor/logs"
    rm -fr "${PROJECT_DIR}/Plot/log_file"
    rm -fr "${PROJECT_DIR}/Plot/plots/variables_dataVmc/"*_part_*.root
    rm -fr "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/Checks/Bias_nominal/condor_logs"
    rm -fr "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/impact_logs"
else
    echo_skip "HiggsZaAna: clean previous outputs"
fi

# The story after HDNA parquets have been produced.

if [[ "$RUN_HIGGSZA_MERGE_PARQUET" == "1" ]]; then
    echo_step "HiggsZaAna: merge parquet"

    activate_conda_env "hza_ana" \
        "$MINICONDA_SETUP" \
        "$MINICONDA_CONDA_SH" \
        "$ANACONDA_SETUP" \
        "$ANACONDA_CONDA_SH" \
        "${CONDA_EXE:+${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh}" \
        "${HOME}/miniconda3/etc/profile.d/conda.sh" \
        "${HOME}/anaconda3/etc/profile.d/conda.sh" \
        "/usr/etc/profile.d/conda.sh"

    # It only works when running the 3_merge_parquet.sh on HiggsDNA directory, otherwise the relative paths in the script will be messed up. So we cd to HiggsDNA directory before running the script, and then cd back to project directory after the script is done.
    cd "${PROJECT_DIR}/HiggsDNA"
    bash "scripts/3_merge_parquet.sh"
else
    echo_skip "HiggsZaAna: merge parquet"
fi

# p2root
if [[ "$RUN_HIGGSZA_P2ROOT" == "1" ]]; then
    echo_step "HiggsZaAna: p2root"

    activate_conda_env "higgs-alp-ana" \
        "$ANACONDA_SETUP" \
        "$ANACONDA_CONDA_SH" \
        "$MINICONDA_SETUP" \
        "$MINICONDA_CONDA_SH" \
        "${CONDA_EXE:+${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh}" \
        "${HOME}/anaconda3/etc/profile.d/conda.sh" \
        "${HOME}/miniconda3/etc/profile.d/conda.sh" \
        "/usr/etc/profile.d/conda.sh"

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
    deactivate_conda_env_if_active
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
    echo_step "AN updated and pushed"
else
    echo_skip "Update AN"
fi

cd "$PROJECT_DIR"
