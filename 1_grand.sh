#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"

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
cd "${PROJECT_DIR}/Parquet2Rootfile"
bash 1_run_P2Root.sh
bash 2_prepare_rootfile.sh

# Train MVA
cd "${PROJECT_DIR}/HZaMVA/scripts"
python3 1_make_sideband_reweight.py
bash 2_train.sh
bash 3_save_model.sh

# p2root for MVA Score
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
cd "${PROJECT_DIR}/Plot"
bash 1_runPlot.sh
