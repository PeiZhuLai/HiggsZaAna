#!/usr/bin/env bash
# Lepton-cleaned (dR>0.3) re-production of the HZa_merged friend parquet, CERN part.
#
# Several copies of this worker share one tag queue (tags_cern.txt). A tag is
# claimed with mkdir, which is atomic on AFS, so two workers never run the same
# tag. Per tag:
#   1. run the existing driver (run_merged_{data,bkg}_ml_friend.sh) for that one
#      tag into the NEW output base -- it blocks until HiggsDNA has finished or
#      retired every job;
#   2. merge by scanning the disk and reconcile, fail-closed (remerge_check.py);
#   3. if jobs are missing (retired after 5 attempts), run once more on the SAME
#      analysis_manager.pkl with --unretire_jobs. Rebuilding the pkl instead
#      could renumber the jobs and pair old chunks with different input files;
#   4. append the verdict to status.txt and, on FAIL, to the shared registry.
#
# Why a queue instead of one driver per fixed tag list: run_analysis blocks on a
# tag's slowest job, so a fixed split leaves workers idle behind one long tag.
#
# Usage (normally started by launch.sh as a systemd --user unit):
#   bash queue_worker.sh <worker_id>
set -uo pipefail

WID="${1:?worker id}"
HERE=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/prod_lepclean
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
LOGD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/prod_lepclean
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/MergedAna/logs/PENDING_RESULTS.txt
QUEUE="${QUEUE:-${HERE}/tags_cern.txt}"
CLAIMS="${CLAIMS:-${LOGD}/claims}"
STATUS="${LOGD}/status.txt"

DATA_BASE=/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_lepclean
BKG_BASE="${BKG_BASE:-/eos/cms/store/group/phys_susy/pelai/HZa_merged/parquet_friend_ML_lepclean}"
# bkg base config; the covered-files rerun points it at za_merged_bkgmc_run3_covered.json
BKG_BASE_CONFIG="${BKG_BASE_CONFIG:-metadata/za_merged_bkgmc_run3.json}"

# A systemd unit does not inherit the login environment: set everything here.
ENVROOT=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export PATH="${ENVROOT}/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export CONDA_PREFIX="${ENVROOT}"
export X509_USER_PROXY=/tmp/x509up_u175325
export EOS_MGM_URL=root://eoscms.cern.ch
export HZA_MLPHOTON_STRICT=1
export HZA_STAGE_INPUTS=1
export FPO_DATA="${FPO_DATA:-1}" FPO_BKG="${FPO_BKG:-4}"
export CONDOR_REQ_MEMORY=3000

mkdir -p "${CLAIMS}"
say() { echo "$(date '+%F %T') w${WID} $*" >> "${STATUS}"; }
say "worker start host=$(hostname -s) python=$(command -v python)"

reconcile() {   # kind sample_dir -> prints the verdict line, returns its rc
    local kind="$1" sdir="$2" flag=""
    [ "${kind}" = "data" ] && flag="--data"
    python "${HERE}/remerge_check.py" "${sdir}" ${flag} 2>&1 | tail -1
    return "${PIPESTATUS[0]}"
}

while read -r kind tag; do
    [ -z "${kind}" ] || [ "${kind:0:1}" = "#" ] && continue
    mkdir "${CLAIMS}/${tag}" 2>/dev/null || continue
    log="${LOGD}/${tag}.log"
    if [ "${kind}" = "data" ]; then
        base="${DATA_BASE}"; sdir="${base}/${tag}/${tag}_2024"
        driver=scripts/run_merged_data_ml_friend.sh; fpo="${FPO_DATA}"
        cfgbase=metadata/za_merged_data_ML_2024.json
    else
        base="${BKG_BASE}"; s="${tag#Bkg_}"; s="${s%_2024}"
        sdir="${base}/${tag}/${s}_2024"
        driver=scripts/run_merged_bkg_ml_friend.sh; fpo="${FPO_BKG}"
        cfgbase="${BKG_BASE_CONFIG}"
    fi
    say "START ${kind} ${tag} fpo=${fpo}"
    t0=$(date +%s)
    ( cd "${REPO}" && TAGS="${tag}" OUTDIR_BASE="${base}" FPO="${fpo}" SKIP_MERGED=0 BASE_CONFIG="${cfgbase}" \
        BATCH_SYSTEM=condor bash "${driver}" ) >> "${log}" 2>&1
    rc=$?
    v=$(reconcile "${kind}" "${sdir}"); vrc=$?
    say "RUN1 ${tag} driver_rc=${rc} $(( ($(date +%s)-t0)/60 ))min :: ${v}"

    if [ "${vrc}" -ne 0 ] && echo "${v}" | grep -q 'missing=[1-9]'; then
        cfg="metadata/friend_tmp_configs/${tag}.json"
        ( cd "${REPO}" && export PYTHONPATH="${REPO}" \
            HIGGSDNA_CONDOR_REQ_MEMORY=3000 \
            HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz \
          && python scripts/run_analysis.py --config "${cfg}" --log-level INFO \
               --output_dir "${base}/${tag}" --batch_system condor --merge_outputs \
               --fpo "${fpo}" --unretire_jobs ) >> "${log}" 2>&1
        v=$(reconcile "${kind}" "${sdir}"); vrc=$?
        say "RUN2 ${tag} (unretire) :: ${v}"
    fi

    if [ "${vrc}" -eq 0 ]; then
        say "DONE ${tag}"
    else
        say "FAIL ${tag}"
        echo "$(date '+%F %H:%M') FAIL  prod_lepclean ${tag} :: ${v}" >> "${REG}"
    fi
done < "${QUEUE}"

say "WORKER_DONE"
