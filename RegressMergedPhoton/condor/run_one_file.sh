#!/bin/bash
# Condor worker for MLNanoAOD production on a single MiniAOD file.
#
# Args:
#   $1  input dataset DAS path or root:// url
#   $2  output file basename (no path)
#   $3  output EOS directory
#
# Environment expected from condor:
#   X509_USER_PROXY  -> path to user grid proxy

set -e

INPUT_URL="$1"
OUTPUT_NAME="$2"
OUTPUT_DIR="$3"

echo "==== MLNanoAOD worker ===="
echo "host=$(hostname)"
echo "date=$(date)"
echo "input=${INPUT_URL}"
echo "output=${OUTPUT_DIR}/${OUTPUT_NAME}"
echo "proxy=${X509_USER_PROXY}"
echo "=========================="

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12

cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src
eval `scramv1 runtime -sh`

# Workspace on local scratch
WORKDIR="${TMPDIR:-/tmp/$USER}/mlnano_$$"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"
echo "workdir=${WORKDIR}"

CFG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/test/Prod_MLNanoAOD_run3_mc.py

cmsRun "${CFG}" \
    inputFiles="${INPUT_URL}" \
    outputFile="${OUTPUT_NAME}" \
    maxEvents=-1

# cmsRun output filename: <stem>_numEventN.root when maxEvents>0; without when -1 the original filename is kept
LOCAL_OUT="${OUTPUT_NAME}"
if [ ! -f "${LOCAL_OUT}" ]; then
    # fallback: pick the produced .root
    LOCAL_OUT=$(ls -t *.root | head -1)
fi
echo "local_output=${LOCAL_OUT}"
echo "size=$(stat -c%s "${LOCAL_OUT}" 2>/dev/null)"

mkdir -p "${OUTPUT_DIR}"
cp "${LOCAL_OUT}" "${OUTPUT_DIR}/${OUTPUT_NAME}"
echo "==== done copy to ${OUTPUT_DIR}/${OUTPUT_NAME} ===="

cd /
rm -rf "${WORKDIR}"
