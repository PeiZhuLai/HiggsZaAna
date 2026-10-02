#!/usr/bin/env bash
#
# 2024 DY+jets (DYJetsTo2E/2Mu/2Tau) HiggsDNA re-production WITH the MC overlap veto
# (mc_overlap_tagger.py now matches DYto2E/DYto2Mu/DYto2Tau; before 2026-09-27 no veto was
# applied to these samples). Mirrors the FSR-fix Bkg_MC production exactly (hza_ana env,
# B-mode pack, copy-then-read staging, fpo=1, RequestMemory 3000, za_bkgmc_run3.json
# tagger/systematics) and the v3v8_extrabkg_20260927 driver; only the sample block changes.
#
# MODE=test -> 3 files/sample   -> .../parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyvetoTEST
# MODE=full -> full datasets    -> .../parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024
#
# WHY the last directory component is not "Bkg_MC": CondorManager writes the condor logs to
#   HiggsDNA/eos_logs/<last component of output_dir>; "Bkg_MC" would collide with the
#   validated FSR-fix production's eos_logs/Bkg_MC.
# WHY a separate catalog (metadata/samples/za_dyveto_2024*.json): SampleManager rewrites
#   <catalog>_sample_manager_full.json; zgamma_tutorial_sample_manager_full.json belongs to
#   the validated production.
# WHY HDNA_REDIRECTOR=cms-xrd-global: the INFN redirector failed for some 2024 files
#   ("Unable to locate ... permission denied") in the v3v8 study; the global one worked.
#
# Launched from a systemd --user unit (lxplus el9 kills setsid/nohup on logout).
set -uo pipefail
MODE="${MODE:-test}"
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/dy2024_reprod_20260927
BASE=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
case "$MODE" in
  test) OUT=$BASE/Bkg_MC_dyvetoTEST; CFG=metadata/za_bkgmc_dyveto2024_test.json ;;
  full) OUT=$BASE/Bkg_MC_dyveto2024; CFG=metadata/za_bkgmc_dyveto2024.json ;;
  *) echo "bad MODE=$MODE"; exit 2 ;;
esac

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1   # dasgoclient, voms-proxy-info
set -u
export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$REPO"
export X509_USER_PROXY=/tmp/x509up_u175325
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
export HZA_STAGE_INPUTS=1
export HDNA_REDIRECTOR=root://cms-xrd-global.cern.ch
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
# RESUME=1 keeps analysis_manager.pkl (jobs in the queue are picked up by cluster id).
# The driver re-saves its pickle on /eos/project every loop and a transient FUSE error there
# killed the v3v8 test driver once, so a crash is retried in resume mode.
clean=1; [ "${RESUME:-0}" = "1" ] && clean=0
rc=1
for attempt in 1 2 3 4 5 6; do
  OUTDIR="$OUT" CONFIG="$CFG" FPO=1 CLEAN_ANALYSIS_STATE=$clean UNRETIRE_JOBS=1 \
    CONDOR_REQ_MEMORY=3000 SAMPLE_LIST="DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau" YEARS=2024 DRY_RUN=0 \
    bash scripts/run_ana_bkgmc.sh
  rc=$?
  echo "[$(date '+%F %T')] attempt $attempt (clean=$clean): run_ana_bkgmc.sh rc=$rc (NOT a verdict; check summaries/chunks)"
  [ $rc -eq 0 ] && break
  [ -s "$OUT/analysis_manager.pkl" ] || { echo "no pickle to resume from -- giving up"; break; }
  clean=0
  sleep 120
done
echo "DRIVER_DONE rc=$rc $(date '+%F %T')" > "$D/logs/driver_${MODE}.done"
exit $rc
