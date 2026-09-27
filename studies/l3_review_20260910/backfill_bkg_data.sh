#!/usr/bin/env bash
#
# Backfill the bkg + data re-production after the 2026-09-12 proxy-expiry wipeout.
#
# WHAT HAPPENED
#   Both stages were submitted 20:40-21:16 on 2026-09-12 while the grid proxy had
#   ~1.4 h left. The proxy was renewed at 21:03, but the jobs had already been
#   given the OLD one, and a proxy is shipped with the job -- renewing on the
#   submit node does nothing for jobs that already hold a copy. Measured outcome:
#       started 20:xx-21:xx : 6990 jobs, all at once
#       exit 0              :  899  (383 bkg + 499 data chunks on disk)
#       exit 1              : 2059  xrootd "[3010] Operation not permitted"
#       walltime-removed    : 4515  hung on xrootd until the longlunch 2 h wall
#   The cliff is exact: chunks stopped appearing at 22:11, the old proxy died at
#   ~22:02-22:25, and the first [3010] stack traces are stamped 22:47-22:49. Jobs
#   that finished inside ~70 min succeeded; everything still reading after the
#   proxy died lost authorization mid-file and either threw or hung.
#
#   Nothing is wrong with the code, the samples, or B-mode. B-mode worked -- no
#   ENV_FETCH_FAIL and no EOS I/O storm, which is what it was introduced for.
#
# WHY IT IS A BACKFILL, NOT A RESUBMIT
#   CLEAN_ANALYSIS_STATE=0 keeps analysis_manager.pkl, so the 899 jobs that
#   already produced a chunk are skipped and only the missing ones run. The
#   pickle was built with exactly the sample list and fpo we want, so the usual
#   "stale pickle overrides the CLI" trap works FOR us here.
#   UNRETIRE_JOBS=1 un-retires jobs HiggsDNA gave up on after 5 failures --
#   without it those stay retired and their chunks are lost for good.
#
# RECONFIGURE_JOBS=1 (2026-09-14)
#   The submit files were written on 2026-09-12 with +JobFlavour = "longlunch"
#   (2 h). Overnight 5081 of 6000 jobs were removed for wall time while only 919
#   actually completed -- and the ones that completed needed a median of 5600 s
#   with a p90 of 10495 s, i.e. the job-time distribution straddles the 2 h wall.
#   The template is now "workday" (8 h), but walltime is fixed at SUBMIT time and
#   condor_qedit cannot change it on a live job, so the submit files must be
#   rewritten -- that is what RECONFIGURE_JOBS=1 does.
#   Auth is no longer the problem: last night's failures were
#   "File did not vector_read properly: [ERROR] Operation expired" with zero
#   "Auth failed"/[3010], which is a slow remote read, not the expired proxy.
#
# PRECONDITION
#   A long proxy. Jobs inherit the proxy that exists when they start, so this
#   refuses to run under 24 h left.
#
# USAGE
#   bash backfill_bkg_data.sh              # dry run
#   DRY_RUN=0 bash backfill_bkg_data.sh    # run (both stages, in parallel)
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
STAGES="${STAGES:-bkg data}"

REPO="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
STAGING="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix"
LOGDIR="${D}/logs_fsrfix"
ENVDIR="/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana"

export CONDA_PREFIX="${ENVDIR}"
export PATH="${ENVDIR}/bin:${PATH}"
export PYTHONPATH="${REPO}"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"
# 12000: measured on the 2026-09-15 staging validation, MemoryUsage median
# 7325 MB / p90 9766 MB. The submit files on disk were raised separately
# (fix_memory_in_submit_files.sh) because --reconfigure_jobs does not rewrite
# them; this variable only matters for any job configured from scratch.
export HIGGSDNA_CONDOR_REQ_MEMORY="${HIGGSDNA_CONDOR_REQ_MEMORY:-12000}"
# Copy-then-read. Default in code is on; exported so it is visible in the job ad
# and can be flipped off from one place if it ever regresses.
export HZA_STAGE_INPUTS="${HZA_STAGE_INPUTS:-1}"

left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
if [ "${left:-0}" -lt 86400 ]; then
    echo "[ERROR] proxy under 24 h -- this is exactly what wiped out the last attempt." >&2
    echo "        voms-proxy-init --rfc --voms cms -valid 192:00, then rerun." >&2
    exit 1
fi
[ -x "${ENVDIR}/bin/python" ] || { echo "[ERROR] ${ENVDIR}/bin/python missing" >&2; exit 1; }

for st in ${STAGES}; do
    case "$st" in
        bkg)  script=run_ana_bkgmc.sh; outdir="${STAGING}/Bkg_MC" ;;
        data) script=run_ana_data.sh;  outdir="${STAGING}/Data" ;;
        *) echo "unknown stage $st"; exit 1 ;;
    esac
    have=$(find "$outdir" -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l | tr -d ' ')
    echo "--- ${st}: ${have} chunks already on disk, backfilling the rest -> ${LOGDIR}/${st}_backfill.log"
    if [ "$DRY_RUN" = "1" ]; then
        echo "    would run: OUTDIR=${outdir} CLEAN_ANALYSIS_STATE=0 UNRETIRE_JOBS=1 bash ${REPO}/scripts/${script}"
        continue
    fi
    ( cd "${REPO}" && setsid nohup env \
        "CONDA_PREFIX=${ENVDIR}" "PATH=${ENVDIR}/bin:${PATH}" "PYTHONPATH=${REPO}" \
        "HZA_BMODE_PACK=${HZA_BMODE_PACK}" "CONDOR_REQ_MEMORY=${HIGGSDNA_CONDOR_REQ_MEMORY}" \
        "OUTDIR=${outdir}" "CLEAN_ANALYSIS_STATE=0" "UNRETIRE_JOBS=1" "DRY_RUN=0" \
        "RECONFIGURE_JOBS=${RECONFIGURE_JOBS:-0}" "HZA_STAGE_INPUTS=${HZA_STAGE_INPUTS}" \
        bash "scripts/${script}" > "${LOGDIR}/${st}_backfill.log" 2>&1 < /dev/null & )
    echo "    launched"
done
echo "=== done (dry run=${DRY_RUN}) ==="
