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

MASS_TAG="${1:?Usage: $0 <mass_tag> <DAS_dataset>}"
DATASET="${2:?Usage: $0 <mass_tag> <DAS_dataset>}"

OUT_EOS="/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/${MASS_TAG}"
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

# ----- condor args file: input_url output_name output_dir -----
ARGS_TXT="${STAGE}/args.txt"
> "${ARGS_TXT}"
i=0
while IFS= read -r f; do
    i=$((i+1))
    out_name="${MASS_TAG}_$(printf '%04d' ${i}).root"
    echo "root://cms-xrd-global.cern.ch/${f} ${out_name} ${OUT_EOS}" >> "${ARGS_TXT}"
done < "${FILES_TXT}"

# ----- JDL -----
JDL="${STAGE}/job.jdl"
cat > "${JDL}" <<JDLEOF
Universe                = vanilla
Executable              = ${WORKER}
Arguments               = \$(input_url) \$(out_name) \$(out_dir)
should_transfer_files   = NO
use_x509userproxy       = true
x509userproxy           = ${PROXY}
JobBatchName            = "MLPhoton_Regressing_${MASS_TAG}"
Log                     = ${LOG_DIR}/job_\$(Cluster)_\$(Process).log
Output                  = ${LOG_DIR}/job_\$(Cluster)_\$(Process).out
Error                   = ${LOG_DIR}/job_\$(Cluster)_\$(Process).err
+JobFlavour             = "workday"
RequestCpus             = 1
RequestMemory           = 4000
RequestDisk             = 4000000
Queue input_url, out_name, out_dir from ${ARGS_TXT}
JDLEOF

echo "JDL written to: ${JDL}"
echo "Submit with:"
echo "  condor_submit ${JDL}"
