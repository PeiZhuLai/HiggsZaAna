#!/bin/bash
# Submit a single mass point (one MiniAOD dataset) to condor for MLNanoAOD production.
#
# Usage:
#   submit_mass_point.sh <mass_tag> <DAS_dataset_path>
#
# Example:
#   submit_mass_point.sh mA_M0p4 \
#     /HZa-Zto2L-ato2G_Par-M-0p4_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM

set -e

MASS_TAG="${1:?Usage: $0 <mass_tag> <DAS_dataset> [mc|data] [year]}"
DATASET="${2:?Usage: $0 <mass_tag> <DAS_dataset> [mc|data] [year]}"
MODE="${3:-mc}"   # mc (default) or data
YEAR="${4:-2024}" # data era (ignored for mc)
FILES_PER_JOB="${FILES_PER_JOB:-1}"  # how many input files cmsRun consumes per condor job

# Output base overridable via OUT_EOS_BASE. Default = pelai group EOS (NPS/project
# EOS were full; HZa_merged now lives under phys_susy/pelai, matching the friend
# producer's MLBASE so orchestrate_data_ml_production.sh finds these files).
OUT_EOS_BASE="${OUT_EOS_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD}"
OUT_EOS="${OUT_EOS_BASE}/${MASS_TAG}"
STAGE="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor/stage/${MASS_TAG}"
LOG_DIR="${STAGE}/logs"
WORKER="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor/run_one_file.sh"
PROXY="${X509_USER_PROXY_AFS:-/afs/cern.ch/user/p/pelai/.x509up_pelai}"

if [ ! -f "${PROXY}" ]; then
    echo "ERROR: voms proxy not found at ${PROXY}."
    echo "Run: voms-proxy-init -voms cms -valid 192:00 && cp /tmp/x509up_u\$(id -u) ${PROXY} && chmod 600 ${PROXY}"
    exit 1
fi

mkdir -p "${OUT_EOS}" "${STAGE}" "${LOG_DIR}"

echo "Mass tag : ${MASS_TAG}"
echo "Dataset  : ${DATASET}"
echo "EOS out  : ${OUT_EOS}"
echo "Stage    : ${STAGE}"
echo

# ----- file list from DAS -----
source /cvmfs/cms.cern.ch/cmsset_default.sh > /dev/null
FILES_TXT="${STAGE}/files.txt"
if [ ! -s "${FILES_TXT}" ]; then
    dasgoclient -query "file dataset=${DATASET}" > "${FILES_TXT}"
fi
NFILES=$(wc -l < "${FILES_TXT}")
echo "Got ${NFILES} files in ${FILES_TXT}"

# ----- condor args file: <inputs-list-file> <output_name> <output_dir> -----
# Per-group URL lists are written to stage/<tag>/inputs/group_NNNN.txt with one
# URL per line. args.txt has exactly 3 whitespace-separated tokens per line and
# NO commas anywhere — HTCondor treats both spaces and commas as field
# separators, so comma-joined URLs in args.txt break the queue parser.
INPUTS_DIR="${STAGE}/inputs"
mkdir -p "${INPUTS_DIR}"
rm -f "${INPUTS_DIR}"/group_*.txt   # clean stale groups

ARGS_TXT="${STAGE}/args.txt"
> "${ARGS_TXT}"
i=0          # cumulative file counter
gi=0         # group counter (= condor job index)
group_file=""
while IFS= read -r f; do
    i=$((i+1))
    if [ -z "${group_file}" ]; then
        gi=$((gi+1))
        group_file="${INPUTS_DIR}/group_$(printf '%04d' ${gi}).txt"
        > "${group_file}"
    fi
    echo "root://cms-xrd-global.cern.ch/${f}" >> "${group_file}"
    if [ $((i % FILES_PER_JOB)) -eq 0 ]; then
        out_name="${MASS_TAG}_$(printf '%04d' ${gi}).root"
        echo "${group_file} ${out_name} ${OUT_EOS}" >> "${ARGS_TXT}"
        group_file=""
    fi
done < "${FILES_TXT}"
# flush remainder
if [ -n "${group_file}" ]; then
    out_name="${MASS_TAG}_$(printf '%04d' ${gi}).root"
    echo "${group_file} ${out_name} ${OUT_EOS}" >> "${ARGS_TXT}"
fi
echo "Will submit ${gi} condor jobs (${FILES_PER_JOB} files per job)"

# ----- JDL -----
JDL="${STAGE}/job.jdl"
cat > "${JDL}" <<JDLEOF
Universe                = vanilla
Executable              = ${WORKER}
Arguments               = \$(input_url) \$(out_name) \$(out_dir) ${MODE} ${YEAR}
should_transfer_files   = NO
use_x509userproxy       = true
x509userproxy           = ${PROXY}
JobBatchName            = "MLPhoton_Regressing_${MASS_TAG}"
Log                     = ${LOG_DIR}/all_jobs.log
Output                  = ${LOG_DIR}/job_\$(Process).out
Error                   = ${LOG_DIR}/job_\$(Process).err
+JobFlavour             = "workday"
RequestCpus             = 1
RequestMemory           = 4000
RequestDisk             = 4000000
max_materialize         = 40
Queue input_url out_name out_dir from ${ARGS_TXT}
JDLEOF

echo "JDL written to: ${JDL}"
echo "Submit with:"
echo "  condor_submit ${JDL}"
