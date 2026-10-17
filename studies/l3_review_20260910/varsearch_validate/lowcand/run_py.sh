#!/usr/bin/env bash
# Condor payload: run one python command of the varsearch validation with the packed higgs-alp-ana env (xrdcp'd
# tarball, no EOS conda activation; see Parquet2Rootfile/Condor/2_submit_dyrescore.sub for why).
# Usage: run_py.sh <script> [args...]
set -uo pipefail
TAR=root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/higgs-alp-ana.tar.gz
W=${_CONDOR_SCRATCH_DIR:-$PWD}/pyenv; mkdir -p $W/env
echo "[JOB] $* host=$(hostname) $(date '+%F %T')"
for t in 1 2 3; do xrdcp -f -s $TAR $W/env.tar.gz && break; sleep 20; done
[ -s $W/env.tar.gz ] || { echo "[ERROR] tarball fetch failed"; exit 2; }
tar -xzf $W/env.tar.gz -C $W/env && rm -f $W/env.tar.gz
export CONDA_PREFIX=$W/env PATH=$W/env/bin:/usr/bin:/bin LD_LIBRARY_PATH=$W/env/lib OMP_NUM_THREADS=1
[ -x $W/env/bin/conda-unpack ] && $W/env/bin/conda-unpack
$W/env/bin/python3 "$@"; rc=$?
rm -rf $W
echo "[JOB] rc=$rc $(date '+%F %T')"
exit $rc
