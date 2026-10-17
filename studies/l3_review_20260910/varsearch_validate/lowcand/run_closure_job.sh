#!/usr/bin/env bash
# Condor payload for the closure of a low-mass candidate (same scripts as the local runs, AFS code + cvmfs CMSSW):
#   run_closure_job.sh scan <model> <mA> <cut>                 one cut of scan_closure.sh (DROP=1, NPAR=1)
#   run_closure_job.sh boot <model> <mA> <cut> <seed> <card>    one replica of boot_closure.sh
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
cd $V || exit 3
echo "[JOB] $* host=$(hostname) $(date '+%F %T')"
case $1 in
  scan) DROP=1 NPAR=1 bash scan_closure.sh $2 $3 "$4" ;;
  boot) BOOT_SEEDS=$5 bash boot_closure.sh $2 $3 $4 30 $6 1 ;;
  *) echo "mode must be scan or boot"; exit 2 ;;
esac
rc=$?; echo "[JOB] rc=$rc $(date '+%F %T')"; exit $rc
