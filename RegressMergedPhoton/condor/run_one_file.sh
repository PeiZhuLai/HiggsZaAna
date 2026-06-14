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

INPUT_ARG="$1"
OUTPUT_NAME="$2"
OUTPUT_DIR="$3"
MODE="${4:-mc}"   # mc (default) or data
YEAR="${5:-2024}" # only used for data: 2022preEE / 2022postEE / 2023preBPix / 2023postBPix / 2024

# INPUT_ARG is either:
#   (a) a single root://... URL, or
#   (b) a path to a text file containing one root://... URL per line (FILES_PER_JOB>1)
# Resolve to a comma-separated URL list for cmsRun inputFiles=.
if [[ "${INPUT_ARG}" == root://* ]] || [[ "${INPUT_ARG}" == http* ]]; then
    INPUT_URL="${INPUT_ARG}"
elif [ -f "${INPUT_ARG}" ]; then
    INPUT_URL=$(paste -sd, "${INPUT_ARG}")
else
    echo "ERROR: input arg is neither a URL nor a readable file: ${INPUT_ARG}"
    exit 1
fi

echo "==== MLNanoAOD worker ===="
echo "host=$(hostname)"
echo "date=$(date)"
echo "input_arg=${INPUT_ARG}"
echo "n_input_files=$(echo "${INPUT_URL}" | tr ',' '\n' | wc -l)"
echo "output=${OUTPUT_DIR}/${OUTPUT_NAME}"
echo "mode=${MODE} year=${YEAR}"
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

if [ "${MODE}" = "data" ]; then
    CFG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/test/Prod_MLNanoAOD_run3_data.py
else
    CFG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/test/Prod_MLNanoAOD_run3_mc.py
fi

# Run cmsRun ONCE PER INPUT URL — one output per input file.
# Each output is named "<tag-prefix>_<MiniAOD-UUID>.root" so downstream HiggsDNA
# can pair it 1:1 with the NanoAODv15 child by querying DAS for the same UUID.
url_list_csv="${INPUT_URL}"
tag_prefix=$(basename "${OUTPUT_NAME}" .root | sed -E 's/_[0-9]+$//')   # strip trailing _NNNN
# Resilience: transient XrdAdaptor file-open failures via the global redirector
# are common at scale. Retry each file a few times; if it still fails, SKIP it
# (do NOT abort the whole job — that would throw away the other good outputs).
# Missing UUIDs are recovered by a resubmit pass (output-presence is the ledger).
i=0; nok=0; nfail=0; FAILED_URLS=""
for url in ${url_list_csv//,/ }; do
    i=$((i+1))
    # MiniAOD basename without extension is its UUID-like identifier
    uuid=$(basename "${url}" .root)
    sub_out="${tag_prefix}_${uuid}.root"
    echo "---- cmsRun $i: ${url} -> ${sub_out}"
    ok=0
    for attempt in 1 2 3; do
        if cmsRun "${CFG}" \
            inputFiles="${url}" \
            outputFile="${sub_out}" \
            maxEvents=-1 \
            year="${YEAR}"; then
            ok=1; break
        fi
        echo "WARN: cmsRun attempt ${attempt}/3 failed for ${uuid}; cleaning partial + backoff"
        rm -f "${sub_out}"
        sleep $((attempt * 30))
    done
    if [ "${ok}" -eq 1 ]; then
        nok=$((nok+1))
    else
        echo "ERROR: cmsRun failed after 3 attempts for ${url} (uuid ${uuid}) — skipping"
        nfail=$((nfail+1)); FAILED_URLS="${FAILED_URLS} ${url}"
    fi
done
echo "cmsRun summary: ok=${nok} fail=${nfail} of ${i}"
[ -n "${FAILED_URLS}" ] && echo "FAILED_URLS:${FAILED_URLS}"
# A job that produced nothing is a real failure → exit 1 so condor marks it.
if [ "${nok}" -eq 0 ]; then
    echo "ERROR: no outputs produced by this job (all ${i} inputs failed)"
    exit 1
fi

# Copy each per-URL output to EOS, preserving the 1:1 UUID mapping.
# On batch nodes the EOS home/user fuse mount (eoshome-p) may not be writable,
# so fall back to xrdcp via the matching MGM if a plain cp fails.
EOS_MGM=""
case "${OUTPUT_DIR}" in
    /eos/home-*) inst=$(echo "${OUTPUT_DIR}" | sed -E 's#^/eos/(home-[a-z]).*#\1#'); EOS_MGM="root://eos${inst}.cern.ch" ;;
    /eos/user/*) EOS_MGM="root://eosuser.cern.ch" ;;
    /eos/project*) EOS_MGM="root://eosproject.cern.ch" ;;
esac

mkdir -p "${OUTPUT_DIR}" 2>/dev/null
if [ -n "${EOS_MGM}" ]; then
    xrdfs "${EOS_MGM}" mkdir -p "${OUTPUT_DIR}" 2>/dev/null || true
fi

for sub_root in "${tag_prefix}"_*.root; do
    [ -f "${sub_root}" ] || continue
    dst="${OUTPUT_DIR}/${sub_root}"
    if cp "${sub_root}" "${dst}" 2>/dev/null; then
        echo "copied (cp): ${dst} ($(stat -c%s "${sub_root}") bytes)"
    elif [ -n "${EOS_MGM}" ] && xrdcp -f "${sub_root}" "${EOS_MGM}/${dst}"; then
        echo "copied (xrdcp): ${EOS_MGM}/${dst} ($(stat -c%s "${sub_root}") bytes)"
    else
        echo "ERROR: copy failed (cp and xrdcp): ${sub_root} -> ${dst}"
        exit 1
    fi
done
echo "==== done all UUID outputs to ${OUTPUT_DIR}/ ===="

cd /
rm -rf "${WORKDIR}"
