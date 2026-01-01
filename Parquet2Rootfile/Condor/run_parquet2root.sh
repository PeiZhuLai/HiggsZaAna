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

# -------- enable conda in THIS shell (no subshell) --------
# 临时关闭 nounset，避免 conda 激活脚本里的未定义变量触发错误
nounset_was_on=0
case $- in *u*) nounset_was_on=1; set +u ;; esac

# 让 'conda' 成为当前 shell 的函数
export PATH="/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH"
if [[ -f /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh ]]; then
  # shellcheck disable=SC1091
  source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
else
  # shellcheck disable=SC1091
  source /eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh || true
fi

# 激活目标环境（失败也不立刻退出，后面再兜底）
conda activate higgs-alp-ana || true

# 激活完成后，按需恢复 nounset
(( nounset_was_on )) && set -u

# 当前 shell 中定位 python
PY_BIN="$(command -v python || true)"

# -------- fallbacks --------
if [[ -z "${PY_BIN:-}" || ! -x "$PY_BIN" ]]; then
  # 1) 直接尝试环境的绝对路径（你的环境路径来自之前报错信息）
  if [[ -x /eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python ]]; then
    PY_BIN="/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python"
  fi
fi

if [[ -z "${PY_BIN:-}" || ! -x "$PY_BIN" ]]; then
  # 2) CVMFS 兜底（可选）
  if [[ -f /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh ]]; then
    # shellcheck disable=SC1091
    source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
    PY_BIN="$(command -v python3 || true)"
  fi
fi

if [[ -z "${PY_BIN:-}" || ! -x "$PY_BIN" ]]; then
  echo "[FATAL] python not found. Tried: conda activate, absolute env path, CVMFS." >&2
  exit 127
fi

# （可选）打印诊断信息，稳定后可注释
"$PY_BIN" -V || true
"$PY_BIN" -c 'import sys,platform; print("PYTHON:",sys.executable); print("PLATFORM:",platform.platform())' || true

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
