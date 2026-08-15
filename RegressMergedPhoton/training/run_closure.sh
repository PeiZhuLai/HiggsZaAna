#!/bin/bash
#
# One-shot MLPhoton closure: dump N events of a signal MiniAOD through the
# online producer + the RecHit dumper, then check that the offline python
# preprocessing reproduces the online MLPhotons.
#
# Input : one MiniAOD file (looked up on DAS if not given)
# Output: <DUMP_DIR>/closure_<TAG>_<N>ev.root      (ROOT, EOS)
#         <DUMP_DIR>/closure_<TAG>_<N>ev.json      (summary)
#         log on stdout
#
# Usage:
#     bash run_closure.sh [TAG] [N_EVENTS] [MINIAOD_FILE]
#     bash run_closure.sh M1 200
#
# Needs a valid grid proxy for the DAS lookup / remote read:
#     voms-proxy-init --rfc --voms cms -valid 192:00
#
set -o pipefail

TAG=${1:-M1}
N_EVENTS=${2:-200}
INFILE=${3:-}

PROJECT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
CMSSW=${PROJECT}/CMSSW_15_0_14
CFG=${CMSSW}/src/RegressMergedPhoton/RecHitDumper/test/closure_cfg.py
TRAINDIR=${PROJECT}/RegressMergedPhoton/training
DUMP_DIR=/eos/home-p/pelai/HZa/root_rechit
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hzgdna/bin/python   # onnxruntime + uproot

OUT_ROOT=${DUMP_DIR}/closure_${TAG}_${N_EVENTS}ev.root
OUT_JSON=${DUMP_DIR}/closure_${TAG}_${N_EVENTS}ev.json

# TAG "M0p5" -> dataset mass token "0p5"
DATASET="/HZa-Zto2L-ato2G_Par-M-${TAG#M}_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM"

export X509_USER_PROXY="/tmp/x509up_u$(id -u)"

mkdir -p "${DUMP_DIR}"

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cd "${CMSSW}/src" || exit 1
eval $(scramv1 runtime -sh) || exit 1
set -u

if [ -z "${INFILE}" ]; then
    echo "[closure] DAS lookup: ${DATASET}"
    INFILE=$(dasgoclient -query "file dataset=${DATASET}" | head -1)
    if [ -z "${INFILE}" ]; then
        echo "[closure] FATAL: no files found on DAS for ${DATASET}"
        exit 1
    fi
fi
echo "[closure] input : ${INFILE}"
echo "[closure] output: ${OUT_ROOT}"

# outFile (not outputFile): VarParsing rewrites outputFile to <stem>_numEvent<N>.root
cmsRun "${CFG}" \
    inputFiles="${INFILE}" \
    outFile="${OUT_ROOT}" \
    maxEvents="${N_EVENTS}" || exit 1

echo "[closure] comparing offline vs online ..."
"${PY}" "${TRAINDIR}/closure_test.py" \
    --input "${OUT_ROOT}" \
    --max-events "${N_EVENTS}" \
    --json "${OUT_JSON}"
