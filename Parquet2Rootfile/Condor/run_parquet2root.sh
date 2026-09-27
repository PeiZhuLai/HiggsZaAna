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

# -------- conda environment setup, matching Plot/Condor behavior --------
PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PY_BIN="${PY_BIN:-python3}"
SETUP_CONDA_ENV="${SETUP_CONDA_ENV:-auto}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
LOCALIZE_CONDA_ENV="${LOCALIZE_CONDA_ENV:-auto}"
LOCAL_CONDA_DIR="${LOCAL_CONDA_DIR:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/higgs-alp-ana-conda}"
# 2026-09-22: the packed env lives on EOS, not in Plot/Condor/env_cache -- the packer
# deletes the local copy after uploading. Defaulting CONDA_TARBALL to that AFS path made
# every job exit 2 before running anything ("tarball is missing"). Empty now means
# "xrdcp it from EOS", exactly like Plot/Condor/run_dataVmc_condor_job.sh does; set
# CONDA_TARBALL explicitly only for an offline/debug local copy.
CONDA_TARBALL_URL="${CONDA_TARBALL_URL:-root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/higgs-alp-ana.tar.gz}"
CONDA_TARBALL="${CONDA_TARBALL:-}"

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

  rm -rf "$LOCAL_CONDA_DIR"
  mkdir -p "$LOCAL_CONDA_DIR"

  local _pack
  if [[ -n "$CONDA_TARBALL" ]]; then
    if [[ ! -s "$CONDA_TARBALL" ]]; then
      echo "[ERROR] CONDA_TARBALL set but missing: $CONDA_TARBALL" >&2
      exit 2
    fi
    _pack="$CONDA_TARBALL"
    echo "[ENV] Extract conda tarball from $_pack to $LOCAL_CONDA_DIR"
  else
    _pack="${LOCAL_CONDA_DIR}.tar.gz"
    echo "[ENV] Fetch conda tarball: $CONDA_TARBALL_URL"
    local _t
    for _t in 1 2 3; do
      xrdcp -f -s "$CONDA_TARBALL_URL" "$_pack" && break
      echo "[ENV] xrdcp attempt $_t failed, retrying in 15s" >&2
      sleep 15
    done
    if [[ ! -s "$_pack" ]]; then
      echo "[ERROR] Could not fetch conda tarball: $CONDA_TARBALL_URL" >&2
      echo "[ERROR] Run Plot/Condor/pack_conda_env_for_condor.sh and upload it, or set CONDA_TARBALL to a local copy." >&2
      exit 2
    fi
  fi

  tar -xzf "$_pack" -C "$LOCAL_CONDA_DIR"
  # only remove the copy we fetched ourselves, never a user-supplied local file
  [[ -z "$CONDA_TARBALL" ]] && rm -f "$_pack"

  # 2026-09-22: export PATH BEFORE conda-unpack, not after. conda-pack writes
  # conda-unpack with a `#!/usr/bin/env python` shebang, so it needs a `python` on PATH
  # to run at all. This used to work by accident: with SETUP_CONDA_ENV=auto the job had
  # already activated the EOS env, which put a python on PATH. With the EOS activation
  # removed there is none, and conda-unpack died with
  #     /usr/bin/env: 'python': No such file or directory
  # taking the whole job down with return value 127 under `set -e`.
  export CONDA_PREFIX="$LOCAL_CONDA_DIR"
  export PATH="$LOCAL_CONDA_DIR/bin:$PATH"
  export LD_LIBRARY_PATH="$LOCAL_CONDA_DIR/lib:${LD_LIBRARY_PATH:-}"
  PY_BIN="$LOCAL_CONDA_DIR/bin/python3"

  if [[ -x "$LOCAL_CONDA_DIR/bin/conda-unpack" ]]; then
    echo "[ENV] Run conda-unpack"
    "$LOCAL_CONDA_DIR/bin/conda-unpack"
  fi
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
