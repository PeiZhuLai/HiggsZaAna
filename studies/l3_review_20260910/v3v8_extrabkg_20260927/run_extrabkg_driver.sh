#!/usr/bin/env bash
#
# V3 (ttbar, ttbar+gamma, ttbar+gammagamma) and V8 (Z+gammagamma) reviewer study:
# HiggsDNA preselection of the extra 2024 background samples, mirroring the FSR-fix
# Bkg_MC production (studies/l3_review_20260910/submit_fpo1_production.sh) exactly:
#   hza_ana env, B-mode pack, copy-then-read staging, fpo=1, RequestMemory 3000,
#   CLEAN_ANALYSIS_STATE=1, same tagger config (za_bkgmc_run3.json) with only the
#   sample block changed.
#
# MODE=test  -> 2 files per sample, output .../parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkgTEST
# MODE=full  -> full TTG / TTGG / ZGG + 10% TTto2L2Nu subset,
#               output .../parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkg2024
#
# WHY the last directory component is not "Bkg_MC"
#   CondorManager puts the condor logs in HiggsDNA/eos_logs/<last component of output_dir>,
#   so ".../Bkg_MC" would write into the existing eos_logs/Bkg_MC of the FSR-fix production.
#
# WHY a separate catalog file (metadata/samples/za_extrabkg_2024*.json)
#   SampleManager.update_catalog() rewrites <catalog>_sample_manager_full.json. Running
#   with zgamma_tutorial.json would overwrite zgamma_tutorial_sample_manager_full.json,
#   which belongs to the validated production. The same four entries are also appended
#   to zgamma_tutorial.json (backup .bak_ttzgg_20260927) for the record.
#
# Launched as a systemd --user unit (lxplus el9 kills setsid/nohup on logout).
set -uo pipefail
MODE="${MODE:-test}"
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927
BASE=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_extraBkg
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
case "$MODE" in
  test) OUT=$BASE/Bkg_MC_extraBkgTEST; CFG=metadata/za_bkgmc_extrabkg2024_test.json ;;
  full) OUT=$BASE/Bkg_MC_extraBkg2024; CFG=metadata/za_bkgmc_extrabkg2024.json ;;
  *) echo "bad MODE=$MODE"; exit 2 ;;
esac

# set -u must be off while sourcing: cmsset_default.sh dereferences unset variables and,
# in the non-interactive systemd environment, aborts the script with an empty log.
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1   # dasgoclient, voms-proxy-info
set -u
export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$REPO"
export X509_USER_PROXY=/tmp/x509up_u175325
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
export HZA_STAGE_INPUTS=1
export HIGGSDNA_CONDOR_REQ_MEMORY=3000

echo "[$(date '+%F %T')] MODE=$MODE OUT=$OUT CFG=$CFG host=$(hostname)"
left=$(voms-proxy-info -file "$X509_USER_PROXY" -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
[ "${left:-0}" -lt 86400 ] && { echo "[ERROR] proxy under 24 h"; exit 1; }
# The repo copy is what ships with each job (jobs.py GRID_PROXY); /tmp being fresh is not enough.
cp -p "$X509_USER_PROXY" "$REPO/x509up_u175325"
openssl x509 -noout -enddate -in "$REPO/x509up_u175325"

mkdir -p "$OUT"
cd "$REPO"
# RESUME=1: keep analysis_manager.pkl (jobs already in the queue are picked up by
# cluster id). The first attempt of a fresh run uses CLEAN_ANALYSIS_STATE=1.
#
# 2026-09-27 14:38: the test driver died on a transient /eos/project FUSE error while
# re-writing analysis_manager_temp.pkl ("mv: cannot stat ... _temp.pkl", then
# FileNotFoundError on open), with all 8 jobs still idle in the queue and the files
# present a moment later. The driver re-saves the pickle on EOS every loop, so over a
# multi-hour run this will recur; restart it in resume mode instead of losing the run.
clean=1; [ "${RESUME:-0}" = "1" ] && clean=0
rc=1
for attempt in 1 2 3 4 5 6; do
  OUTDIR="$OUT" CONFIG="$CFG" FPO=1 CLEAN_ANALYSIS_STATE=$clean UNRETIRE_JOBS=1 \
    CONDOR_REQ_MEMORY=3000 SAMPLE_LIST="TTto2L2Nu,TTG,TTGG_Run3,ZGG" YEARS=2024 DRY_RUN=0 \
    bash scripts/run_ana_bkgmc.sh
  rc=$?
  echo "[$(date '+%F %T')] attempt $attempt (clean=$clean): run_ana_bkgmc.sh rc=$rc (NOT a verdict; check summaries/chunks)"
  [ $rc -eq 0 ] && break
  [ -s "$OUT/analysis_manager.pkl" ] || { echo "no pickle to resume from -- giving up"; break; }
  clean=0
  sleep 120
done
echo "DRIVER_DONE rc=$rc" > "$D/logs/driver_${MODE}.done"
exit $rc
