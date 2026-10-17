#!/usr/bin/env bash
# Condor job: re-score every era present of one sample dir (data / DY too) with a given low-mass model (rescore_selected.py --low-model; copy of run_rescore_sig.sh;
# low model). Inputs are xrdcp'd to the job scratch first (FUSE reads made a local test I/O bound: 4.4 min wall,
# 30 s CPU); the output is xrdcp'd back to <out>.partial with an adler32 end-to-end check, then renamed. Eras whose
# output already exists are skipped. Conda env from the packed tarball (LOCALIZE=1, no EOS activation, see the
# notes in Parquet2Rootfile/Condor/2_submit_dyrescore.sub).
# Usage: run_rescore_any.sh <sample> <eos out base> <low model name>   exit 0 only when every present era has its output
set -uo pipefail
SAMPLE=$1; OUT=$2; LOWM=$3; HIGHM=${4:-high}
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
SRC=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix
XR=root://eosuser.cern.ch/
# xrootd namespace of /eos/home-p/pelai is /eos/user/p/pelai (the FUSE alias does not exist on the server)
xp() { echo "$1" | sed 's#^/eos/home-p/pelai/#/eos/user/p/pelai/#'; }
XSRC=$(xp $SRC); XOUT=$(xp $OUT)
TAR=root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/higgs-alp-ana.tar.gz
W=${_CONDOR_SCRATCH_DIR:-$PWD}/rescore_work; mkdir -p $W/in/$SAMPLE $W/out $W/env
echo "[JOB] $SAMPLE host=$(hostname) $(date '+%F %T')"
todo=()
ERAS=$(for t in 1 2 3; do o=$(xrdfs eosuser.cern.ch ls $XSRC/$SAMPLE 2>/dev/null) && { echo "$o"; break; }; sleep 20; done | sed -n 's#.*/\(.*\)\.root$#\1#p')   # single quotes: inside double quotes $# expands
[ -n "$ERAS" ] || { echo "[ERROR] cannot list $XSRC/$SAMPLE"; exit 2; }
echo "[JOB] eras: $ERAS"
for e in $ERAS; do
  sz=$(xrdfs eosuser.cern.ch stat $XOUT/$SAMPLE/$e.root 2>/dev/null | awk '/^Size:/{print $2}')
  [ -n "$sz" ] && [ "$sz" -gt 0 ] && { echo "[SKIP] $e exists"; continue; }; todo+=($e)
done
[ ${#todo[@]} -eq 0 ] && { echo "[JOB] nothing to do"; exit 0; }
echo "[ENV] Fetch conda tarball: $TAR"
for t in 1 2 3; do xrdcp -f -s $TAR $W/env.tar.gz && break; sleep 20; done
[ -s $W/env.tar.gz ] || { echo "[ERROR] tarball fetch failed"; exit 2; }
tar -xzf $W/env.tar.gz -C $W/env && rm -f $W/env.tar.gz
export CONDA_PREFIX=$W/env PATH=$W/env/bin:/usr/bin:/bin LD_LIBRARY_PATH=$W/env/lib
[ -x $W/env/bin/conda-unpack ] && $W/env/bin/conda-unpack
PY=$W/env/bin/python3
rc=0
for e in "${todo[@]}"; do
  ok=0
  for t in 1 2 3; do
    xrdcp -f -s $XR$XSRC/$SAMPLE/$e.root $W/in/$SAMPLE/$e.root || { sleep 20; continue; }
    $PY $V/rescore_selected.py --sample $SAMPLE --era $e --src-base $W/in --out-base $W/out --low-model $LOWM --high-model $HIGHM > $W/py_$e.log 2>&1; prc=$?
    cat $W/py_$e.log
    grep -q "GATE FAIL\|READBACK FAIL" $W/py_$e.log && { echo "[ERROR] $e gate/readback failure, not retried"; break; }
    [ $prc -eq 0 ] || { sleep 20; continue; }
    xrdfs eosuser.cern.ch mkdir -p $XOUT/$SAMPLE
    xrdcp -f -s --cksum adler32:source $W/out/$SAMPLE/$e.root $XR$XOUT/$SAMPLE/$e.root.partial || { sleep 20; continue; }
    xrdfs eosuser.cern.ch mv $XOUT/$SAMPLE/$e.root.partial $XOUT/$SAMPLE/$e.root || { sleep 20; continue; }
    echo "[DONE] $SAMPLE $e adler32=$(xrdadler32 $W/out/$SAMPLE/$e.root | awk '{print $1}')"; ok=1; break
  done
  rm -f $W/in/$SAMPLE/$e.root $W/out/$SAMPLE/$e.root
  [ $ok -eq 1 ] || { echo "[FAILED] $SAMPLE $e"; rc=1; }
done
rm -rf $W
exit $rc
