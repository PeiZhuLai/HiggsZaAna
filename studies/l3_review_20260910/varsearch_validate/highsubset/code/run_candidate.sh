#!/bin/bash
# CERN condor payload for one candidate of the BDT variable search.
#   $1 = index into config.jsonl, $2 = EOS results dir, $3 = model (low|high)
# search_job.py and config.jsonl arrive through transfer_input_files (no AFS read while running).
# The environment comes from cvmfs (LCG_106: numpy 1.26, xgboost 2.0.3), never from the EOS conda.
# Inputs are copied sequentially from EOS into the job scratch dir; the result JSON is copied back to EOS.
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
set -u
IDX=$1; RES=$2; MODEL=$3
SAMPLES=/eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples
WORK=${_CONDOR_SCRATCH_DIR:-$PWD}
mkdir -p $WORK/samples $WORK/out
trap 'rm -rf $WORK/samples $WORK/out' EXIT
if [ "$MODEL" = low ]; then M="1 2 3"; else M="4 5 6 7 8 9 10 15 20 25 30"; fi
files="meta.json data_sb.npz dy.npz train_bkg.npz"
for m in $M; do files="$files test_sig_m$m.npz train_sig_m$m.npz"; done
for f in $files; do
  for t in 1 2 3; do cp $SAMPLES/$f $WORK/samples/ && break; sleep 30; done
  [ -s $WORK/samples/$f ] || { echo "FAILED to copy $f"; exit 3; }
done
echo "[$(date '+%F %T')] inputs copied: $(du -sh $WORK/samples | cut -f1) on $(hostname)"
export OMP_NUM_THREADS=1 NTHREADS=1
python3 search_job.py --samples $WORK/samples --config config.jsonl --index $IDX --out $WORK/out || exit 4
for f in $WORK/out/*.json; do
  b=$(basename $f)
  # copy under a temporary name, then rename: the driver never reads a half-written result
  for t in 1 2 3; do cp $f $RES/.$b.tmp && mv -f $RES/.$b.tmp $RES/$b && break; sleep 30; done
  [ -s $RES/$b ] || { echo "FAILED to copy back $f"; exit 5; }
done
echo "[$(date '+%F %T')] done"
