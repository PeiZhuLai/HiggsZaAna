#!/usr/bin/env bash
#
# Fill the 154 signal jobs that never produced output.
#
# WHY a separate script instead of just rerunning the stage
#   The first driver hung: the CERN schedd bigbird26 stopped accepting queries
#   for ~3 h (SECMAN authorization failures while the schedd sat at 0.75 duty
#   cycle), the driver blocked in do_wait on a condor call that never returned,
#   and it never wrote to its log again after 17:56. Before hanging it had
#   already retired 59 jobs ("submitted 5 times, retiring job"). Left alone it
#   would have exited 0 claiming 100 % complete with 154 chunks missing.
#
# WHAT THIS DOES DIFFERENTLY
#   CLEAN_ANALYSIS_STATE=0  keep the manager state, so finished jobs are not redone
#   UNRETIRE_JOBS=1         pick the retired ones back up
#   Only ~154 jobs are left, so concurrency stays far below the level where the
#   shared conda env on /eos starts returning I/O errors (clean at ~211, 66 %
#   failures at ~1089).
#
# USAGE
#   bash backfill_signal.sh            # dry run
#   DRY_RUN=0 bash backfill_signal.sh  # submit
#
set -uo pipefail

REPO="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
STAGING="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
ENVDIR="/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana"
DRY_RUN="${DRY_RUN:-1}"

export CONDA_PREFIX="${ENVDIR}"
export PATH="${ENVDIR}/bin:${PATH}"
export PYTHONPATH="${REPO}"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"

missing=$(find "${STAGING}" -mindepth 2 -maxdepth 2 -type d -name 'job_*' \
          '!' -exec sh -c 'ls "$1"/*_nominal.parquet >/dev/null 2>&1' _ {} ';' -print 2>/dev/null | wc -l)
echo "job dirs still without a nominal parquet: ${missing}"
left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
echo "dry run: ${DRY_RUN}"

if [ "${DRY_RUN}" = "1" ]; then
    echo "would run: CLEAN_ANALYSIS_STATE=0 UNRETIRE_JOBS=1 OUTDIR=${STAGING} bash ${REPO}/scripts/run_ana_signal.sh"
    exit 0
fi
if [ "${left:-0}" -lt 7200 ]; then
    echo "[ERROR] proxy under 2 h; renew before backfilling." >&2
    exit 1
fi

mkdir -p "${D}/logs_fsrfix"
cd "${REPO}" && env \
    "CONDA_PREFIX=${ENVDIR}" "PATH=${ENVDIR}/bin:${PATH}" "PYTHONPATH=${REPO}" \
    "OUTDIR=${STAGING}" "CLEAN_ANALYSIS_STATE=0" "UNRETIRE_JOBS=1" "RETIRE_JOBS=0" \
    bash scripts/run_ana_signal.sh
echo "backfill driver returned $? -- check artifacts, not this code"
