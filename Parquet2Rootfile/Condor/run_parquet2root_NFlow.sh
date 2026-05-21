#!/bin/bash
# Wrapper for running Parque2Root_BDT_NFlow.py on LXBATCH/HTCondor.
# Usage: run_parquet2root_NFlow.sh <INPUT> <OUTPUT> <CORR> <SPLIT_FLAG>

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

# -------- conda environment setup, matching Plot/Condor behavior --------
PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PY_BIN="${PY_BIN:-python3}"
SETUP_CONDA_ENV="${SETUP_CONDA_ENV:-auto}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
LOCALIZE_CONDA_ENV="${LOCALIZE_CONDA_ENV:-auto}"
LOCAL_CONDA_DIR="${LOCAL_CONDA_DIR:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/higgs-alp-ana-conda}"
CONDA_TARBALL="${CONDA_TARBALL:-${PROJECT_DIR}/Plot/Condor/env_cache/higgs-alp-ana.tar.gz}"

sanitize_python_env_for_conda() {
  if [[ -n "${PYTHONPATH:-}" ]]; then
    echo "[ENV] Unset inherited PYTHONPATH before conda activation"
    unset PYTHONPATH
  fi

  if [[ -n "${PYTHONHOME:-}" ]]; then
    echo "[ENV] Unset inherited PYTHONHOME before conda activation"
    unset PYTHONHOME
  fi
}

activate_conda_env_if_needed() {
  local should_activate=0

  case "$SETUP_CONDA_ENV" in
    1|true|TRUE|yes|YES) should_activate=1 ;;
    0|false|FALSE|no|NO) should_activate=0 ;;
    auto|AUTO)
      if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV_NAME" ]]; then
        should_activate=1
      fi
      ;;
    *)
      echo "[ERROR] SETUP_CONDA_ENV must be auto, 1, or 0; got '$SETUP_CONDA_ENV'" >&2
      exit 2
      ;;
  esac

  if [[ "$should_activate" -ne 1 ]]; then
    return 0
  fi

  echo "[ENV] Activate conda env: $CONDA_ENV_NAME"
  if [[ ! -r "$ANACONDA_SETUP" ]]; then
    echo "[ERROR] Cannot read ANACONDA_SETUP: $ANACONDA_SETUP" >&2
    exit 2
  fi

  echo "[ENV] Source anaconda setup: $ANACONDA_SETUP"
  set +u
  sanitize_python_env_for_conda
  # shellcheck disable=SC1090
  source "$ANACONDA_SETUP"
  conda activate "$CONDA_ENV_NAME"
  set -u

  if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/python3" ]]; then
    echo "[ERROR] Failed to activate conda env: $CONDA_ENV_NAME" >&2
    exit 2
  fi

  PY_BIN="${CONDA_PREFIX}/bin/python3"
  echo "[ENV] active CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
  echo "[ENV] active python=${PY_BIN}"
}

localize_conda_env_if_needed() {
  local source_env="${CONDA_PREFIX:-}"
  local should_localize=0

  case "$LOCALIZE_CONDA_ENV" in
    1|true|TRUE|yes|YES) should_localize=1 ;;
    0|false|FALSE|no|NO) should_localize=0 ;;
    auto|AUTO)
      if [[ -n "$source_env" && "$source_env" == /eos/* ]]; then
        should_localize=1
      fi
      ;;
    *)
      echo "[ERROR] LOCALIZE_CONDA_ENV must be auto, 1, or 0; got '$LOCALIZE_CONDA_ENV'" >&2
      exit 2
      ;;
  esac

  if [[ "$should_localize" -ne 1 ]]; then
    return 0
  fi

  if [[ ! -s "$CONDA_TARBALL" ]]; then
    echo "[ERROR] Conda localization requested, but tarball is missing: $CONDA_TARBALL" >&2
    echo "[ERROR] Run Plot/Condor/pack_conda_env_for_condor.sh once before condor_submit." >&2
    exit 2
  fi

  echo "[ENV] Extract conda tarball from $CONDA_TARBALL to $LOCAL_CONDA_DIR"
  rm -rf "$LOCAL_CONDA_DIR"
  mkdir -p "$LOCAL_CONDA_DIR"
  tar -xzf "$CONDA_TARBALL" -C "$LOCAL_CONDA_DIR"

  if [[ -x "$LOCAL_CONDA_DIR/bin/conda-unpack" ]]; then
    echo "[ENV] Run conda-unpack"
    "$LOCAL_CONDA_DIR/bin/conda-unpack"
  fi

  export CONDA_PREFIX="$LOCAL_CONDA_DIR"
  export PATH="$LOCAL_CONDA_DIR/bin:$PATH"
  export LD_LIBRARY_PATH="$LOCAL_CONDA_DIR/lib:${LD_LIBRARY_PATH:-}"
  PY_BIN="$LOCAL_CONDA_DIR/bin/python3"
}

echo "[ENV] initial PY_BIN=${PY_BIN}"
echo "[ENV] initial CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
echo "[ENV] SETUP_CONDA_ENV=${SETUP_CONDA_ENV}"
echo "[ENV] CONDA_ENV_NAME=${CONDA_ENV_NAME}"
echo "[ENV] ANACONDA_SETUP=${ANACONDA_SETUP}"
echo "[ENV] LOCALIZE_CONDA_ENV=${LOCALIZE_CONDA_ENV}"
echo "[ENV] CONDA_TARBALL=${CONDA_TARBALL}"

activate_conda_env_if_needed
localize_conda_env_if_needed

if [[ "$PY_BIN" != */* ]]; then
  PY_BIN="$(command -v "$PY_BIN" || true)"
fi
if [[ ! -x "$PY_BIN" ]]; then
  echo "[FATAL] python not found after environment setup: $PY_BIN" >&2
  exit 127
fi

"$PY_BIN" -V
"$PY_BIN" -c 'import sys,platform; print("PYTHON:",sys.executable); print("PLATFORM:",platform.platform())'
"$PY_BIN" -c 'import pandas, uproot, ROOT, xgboost; print("IMPORTS: pandas/uproot/ROOT/xgboost OK")'

# -------- Python script path (AFS) --------
PY_SCRIPT="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT_NFlow.py"

# -------- Build command --------
CMD=("$PY_BIN" "$PY_SCRIPT" -i "$INPUT" -o "$OUTPUT")
if [[ "$SPLIT_FLAG" == "1" ]]; then
  CMD+=(--split)
fi
# CORR is kept in the joblist for traceability. The Python script currently infers inputs from -i.

# -------- Retry --------
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
