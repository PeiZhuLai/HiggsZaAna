#!/bin/bash
# Fetch PDF members from the input NanoAOD and compute the per-sample PDF acceptance
# uncertainty for all 70 Run 3 signal samples. One sample at a time, 6 worker processes
# (lxplus foreground limit). Resumable: per-file cache under cache/, finished samples
# have results/<sample>.json. Completion marker: logs/run_all.DONE
set -o pipefail
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/pdf_acceptance_20261002
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
export X509_USER_PROXY=/tmp/x509up_u175325
cd $W
rm -f logs/run_all.DONE
SAMPLES=$(ls -d /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC/mA_M*_20* | xargs -n1 basename | sort -V)
nfail=0
for s in $SAMPLES; do
  if [ -s results/$s.json ] && [ "$1" != "--force" ]; then echo "[skip] $s"; continue; fi
  for try in 1 2 3; do
    $PY pdf_acceptance.py both $s --nproc 6 && break
    echo "[retry $try] $s"; sleep 30
  done
  [ -s results/$s.json ] || { echo "[FAILED] $s"; nfail=$((nfail+1)); }
done
echo "nfail=$nfail" > logs/run_all.DONE
date >> logs/run_all.DONE
