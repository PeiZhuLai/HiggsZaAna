#!/bin/bash
# Wrapper for running Parque2Root_BDT.py on LXBATCH/HTCondor
# Usage: run_parquet2root.sh <INPUT> <OUTPUT> <MA> <CORR> <SPLIT_FLAG>

set -euo pipefail

INPUT="$1"      # /eos/.../merged_*.parquet
OUTPUT="$2"     # /eos/.../*.root
CORR="$3"       # nominal / FNUF_up / ...
SPLIT_FLAG="$4" # "1" for signal (--split), "0" otherwise

# -------- single-thread everything (avoid PSI & oversubscription) --------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export ARROW_NUM_THREADS=1

# -------- enable original conda environment in THIS shell --------
# Keep using the existing Anaconda installation requested for this workflow.
CONDA_BASE="/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3"
CONDA_ENV="higgs-alp-ana"
nounset_was_on=0
case $- in *u*) nounset_was_on=1; set +u ;; esac

export PATH="$CONDA_BASE/bin:$PATH"
if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [[ -f /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh ]]; then
  # shellcheck disable=SC1091
  source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
else
  echo "[FATAL] conda setup script not found under $CONDA_BASE" >&2
  exit 127
fi

conda activate "$CONDA_ENV"

(( nounset_was_on )) && set -u

PY_BIN="$(command -v python || true)"
if [[ -z "${PY_BIN:-}" || ! -x "$PY_BIN" ]]; then
  echo "[FATAL] python not found after activating conda environment: $CONDA_ENV" >&2
  exit 127
fi

"$PY_BIN" -V
"$PY_BIN" -c 'import sys,platform; print("PYTHON:",sys.executable); print("PLATFORM:",platform.platform())'
"$PY_BIN" -c 'import pandas, uproot, ROOT, xgboost; print("IMPORTS: pandas/uproot/ROOT/xgboost OK")'

# -------- Python 脚本路径（AFS）--------
PY_SCRIPT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py"

# -------- 构建命令 --------
CMD=("$PY_BIN" "$PY_SCRIPT" -i "$INPUT" -o "$OUTPUT")
if [[ "$SPLIT_FLAG" == "1" ]]; then
  CMD+=(--split)
fi
# 如需把 CORR 也交给脚本用作内部逻辑：CMD+=(--corr "$CORR")

# -------- 重试机制 --------
MAX_RETRY=3
TRY=1
while (( TRY <= MAX_RETRY )); do
  echo "[$(date)] Attempt ${TRY}: ${CMD[*]}"
  if "${CMD[@]}"; then
    echo "[$(date)] SUCCESS"
    exit 0
  fi
  echo "[$(date)] Failed. Retrying..."
  ((TRY++))
  sleep 10
done

echo "[$(date)] ERROR: Job failed after ${MAX_RETRY} attempts"
exit 1
