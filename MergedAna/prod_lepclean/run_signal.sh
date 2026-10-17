#!/usr/bin/env bash
# Lepton-cleaned signal re-production for HZa_merged, OLD MLPhoton model
# (same model as the data friend, which is what the ROI windows in
# merged_p2root.py assume for MERGED_ML_VERSION=old).
#
# Input : /eos/.../HZa_merged/MLNanoAOD/mA_M0p{1..9}   (old-model friend)
# Output: /eos/.../HZa_merged/parquet_merged_DNA_lepclean_old/Sig_MC_MLNANO_all
# Drives HiggsDNA/scripts/run_merged_signal_v3.sh (its ML_BASE/OUTDIR/GEN_DIR are
# overridden; the "v3" in its name is only the default input), then merges each
# mass point by scanning the disk and reconciles (remerge_check.py).
set -uo pipefail

HERE=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/prod_lepclean
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
STATUS="${LOGD}/status.txt"
OUT=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_merged_DNA_lepclean_old/Sig_MC_MLNANO_all

export PATH="/usr/bin:/bin:/usr/sbin:/sbin"
export X509_USER_PROXY=/tmp/x509up_u175325
export EOS_MGM_URL=root://eoscms.cern.ch
export HZA_MLPHOTON_STRICT=1 HZA_STAGE_INPUTS=1
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
say() { echo "$(date '+%F %T') sig $*" >> "${STATUS}"; }

mkdir -p "${OUT}"
say "START signal -> ${OUT}"
( cd "${REPO}" && ML_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged/MLNanoAOD \
    OUTDIR="${OUT}" GEN_DIR=metadata/_generated/lepclean_old CONDOR_REQ_MEMORY=3000 FPO=1 \
    bash scripts/run_merged_signal_v3.sh ) >> "${LOGD}/signal.log" 2>&1
say "driver rc=$?"

PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
nfail=0
for m in 1 2 3 4 5 6 7 8 9; do
    v=$("${PY}" "${HERE}/remerge_check.py" "${OUT}/mA_MLNANO_M0p${m}_2024" 2>&1 | tail -1)
    rc=$?
    say "M0p${m} :: ${v}"
    echo "${v}" | grep -q '^PASS' || { nfail=$((nfail+1)); echo "$(date '+%F %H:%M') FAIL  prod_lepclean signal M0p${m} :: ${v}" >> "${REG}"; }
done
say "SIGNAL_DONE nfail=${nfail}"
