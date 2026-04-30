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
ENV_MAX_RETRY="${ENV_MAX_RETRY:-6}"
ENV_RETRY_SLEEP="${ENV_RETRY_SLEEP:-60}"
EOS_WARMUP_SLEEP="${EOS_WARMUP_SLEEP:-30}"
nounset_was_on=0
case $- in *u*) nounset_was_on=1; set +u ;; esac

setup_conda_once() {
  export PATH="$CONDA_BASE/bin:$PATH"
  if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
  elif [[ -f /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh ]]; then
    # shellcheck disable=SC1091
    source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh
  else
    echo "[FATAL] conda setup script not found under $CONDA_BASE" >&2
    return 127
  fi

  conda activate "$CONDA_ENV"
}

env_try=1
while (( env_try <= ENV_MAX_RETRY )); do
  echo "[$(date)] Environment setup attempt ${env_try}/${ENV_MAX_RETRY}"
  if setup_conda_once; then
    break
  fi
  if (( env_try == ENV_MAX_RETRY )); then
    echo "[FATAL] failed to activate conda environment after ${ENV_MAX_RETRY} attempts: $CONDA_ENV" >&2
    exit 127
  fi
  echo "[$(date)] Environment setup failed. Waiting ${ENV_RETRY_SLEEP}s for EOS, then retrying..."
  sleep "$ENV_RETRY_SLEEP"
  ((env_try++))
done

(( nounset_was_on )) && set -u

PY_BIN="$(command -v python || true)"
if [[ -z "${PY_BIN:-}" || ! -x "$PY_BIN" ]]; then
  echo "[FATAL] python not found after activating conda environment: $CONDA_ENV" >&2
  exit 127
fi

"$PY_BIN" -V
"$PY_BIN" -c 'import sys,platform; print("PYTHON:",sys.executable); print("PLATFORM:",platform.platform())'
echo "[$(date)] Waiting ${EOS_WARMUP_SLEEP}s before importing EOS-hosted Python packages..."
sleep "$EOS_WARMUP_SLEEP"

import_try=1
while (( import_try <= ENV_MAX_RETRY )); do
  echo "[$(date)] Python import check attempt ${import_try}/${ENV_MAX_RETRY}"
  if "$PY_BIN" -c 'import pandas, uproot, ROOT, xgboost; print("IMPORTS: pandas/uproot/ROOT/xgboost OK")'; then
    break
  fi
  if (( import_try == ENV_MAX_RETRY )); then
    echo "[FATAL] Python import check failed after ${ENV_MAX_RETRY} attempts. EOS may still be unable to read the conda environment files." >&2
    exit 1
  fi
  echo "[$(date)] Import check failed. Waiting ${ENV_RETRY_SLEEP}s for EOS, then retrying..."
  sleep "$ENV_RETRY_SLEEP"
  ((import_try++))
done

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
