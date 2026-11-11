#!/usr/bin/env bash
# Lepton-cleaned signal production for HZa_merged with a chosen MLPhoton model,
# current tagger (lepton cleaning + strict loader) -- the generalization of
# run_signal.sh (which is fixed to the old model and M0p1-M0p9).
#
# Why: the model comparison plots (Plot/scripts/*mlphoton*) compared the old model
# on the current tagger (lepclean_old) with v4 on the August tagger (no lepton
# cleaning), and M1 existed only in the July production. This makes both models
# available on the same tagger, M1 included.
#
#   VERSION=old : MLNanoAOD    -> parquet_merged_DNA_lepclean_old/Sig_MC_MLNANO_all
#   VERSION=v4  : MLNanoAOD_v4 -> parquet_merged_DNA_lepclean_v4cur/Sig_MC_MLNANO_all
#                 (not lepclean_v4: that holds the 09-02 M0p1-only production)
#
# One condor job per nano file (FPO=1), 3000 MB. Each mass point is merged by
# scanning the disk and reconciled fail-closed (remerge_check.py).
#
# Dry-run by default (prints the configs, submits nothing); APPLY=1 to submit.
#   VERSION=v4 bash run_signal_ml.sh M1
#   VERSION=v4 APPLY=1 bash run_signal_ml.sh M1
set -uo pipefail

VERSION="${VERSION:?VERSION=old or v4}"
APPLY="${APPLY:-0}"
MASSES="${*:-M0p1 M0p2 M0p3 M0p4 M0p5 M0p6 M0p7 M0p8 M0p9 M1}"

HERE=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/prod_lepclean
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged
case "${VERSION}" in
    old) ML_BASE="${BASE}/MLNanoAOD";    PROD=parquet_merged_DNA_lepclean_old ;;
    v4)  ML_BASE="${BASE}/MLNanoAOD_v4"; PROD=parquet_merged_DNA_lepclean_v4cur ;;
    *)   echo "unknown VERSION=${VERSION}"; exit 1 ;;
esac
OUT="${BASE}/${PROD}/Sig_MC_MLNANO_all"
STATUS="${LOGD}/status_sig_${VERSION}.txt"

export PATH="/usr/bin:/bin:/usr/sbin:/sbin"
export X509_USER_PROXY=/tmp/x509up_u175325
export EOS_MGM_URL=root://eoscms.cern.ch
export HZA_MLPHOTON_STRICT=1 HZA_STAGE_INPUTS=1
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
say() { echo "$(date '+%F %T') sig_${VERSION} $*" | tee -a "${STATUS}"; }

# ---- pre-flight -----------------------------------------------------------
bad=0
tl=$(voms-proxy-info -file /tmp/x509up_u175325 -timeleft 2>/dev/null || echo 0)
echo "proxy /tmp: $((tl/3600)) h left"; [ "${tl}" -lt 172800 ] && { echo "  !! < 48 h"; bad=1; }
cmp -s /tmp/x509up_u175325 "${REPO}/x509up_u175325" || { echo "  !! repo proxy copy differs from /tmp"; bad=1; }
for m in ${MASSES}; do
    pm="${REPO}/metadata/parent_maps/mA_${m}.json"
    [ -f "${pm}" ] || { echo "  !! no parent map mA_${m}"; bad=1; continue; }
    [ -d "${ML_BASE}/mA_${m}" ] || { echo "  !! no ${ML_BASE}/mA_${m}"; bad=1; continue; }
    [ -e "${OUT}/mA_MLNANO_${m}_2024" ] && { echo "  !! output exists: ${OUT}/mA_MLNANO_${m}_2024"; bad=1; }
    echo "  ${m}: $(python3 -c "import json,sys;print(len(json.load(open(sys.argv[1]))))" "${pm}") nano files (= jobs)"
done
# HiggsDNA batch names are <sample>_<year>, the same for both versions: a mass point
# must not run in both at once (other mass points may)
qn=$(condor_q -af JobBatchName 2>/dev/null)
for m in ${MASSES}; do
    q=$(echo "${qn}" | grep -c "^mA_MLNANO_${m}_2024$" || true)
    [ "${q}" -gt 0 ] && { echo "  !! ${q} mA_MLNANO_${m}_2024 jobs already in the queue (same batch name)"; bad=1; }
done
[ "${bad}" -ne 0 ] && { echo "pre-flight FAILED"; exit 1; }
echo "pre-flight OK  ->  ${OUT}"

if [ "${APPLY}" != "1" ]; then
    ( cd "${REPO}" && DRY_RUN=1 ML_BASE="${ML_BASE}" OUTDIR="${OUT}" \
        GEN_DIR="metadata/_generated/lepclean_${VERSION}" \
        bash scripts/run_merged_signal_v3.sh ${MASSES} )
    echo "(dry-run; APPLY=1 to submit)"
    exit 0
fi

mkdir -p "${OUT}"
say "START ${MASSES} -> ${OUT}"
( cd "${REPO}" && ML_BASE="${ML_BASE}" OUTDIR="${OUT}" \
    GEN_DIR="metadata/_generated/lepclean_${VERSION}" CONDOR_REQ_MEMORY=3000 FPO=1 \
    bash scripts/run_merged_signal_v3.sh ${MASSES} ) >> "${LOGD}/signal_${VERSION}.log" 2>&1
say "driver rc=$?"

PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
nfail=0
for m in ${MASSES}; do
    v=$("${PY}" "${HERE}/remerge_check.py" "${OUT}/mA_MLNANO_${m}_2024" 2>&1 | tail -1)
    say "${m} :: ${v}"
    echo "${v}" | grep -q '^PASS' || { nfail=$((nfail+1)); echo "$(date '+%F %H:%M') FAIL  prod_lepclean signal_${VERSION} ${m} :: ${v}" >> "${REG}"; }
done
say "SIGNAL_DONE nfail=${nfail}"
[ "${nfail}" -eq 0 ] && echo "$(date '+%F %H:%M') DONE  prod_lepclean signal_${VERSION} ${MASSES} -> ${OUT}" >> "${REG}"
