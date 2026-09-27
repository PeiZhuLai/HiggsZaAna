#!/bin/bash
set -euo pipefail

timer_start=$(date +%s)
finish() {
    local status=$?
    local timer_end elapsed hours minutes seconds
    local runtime_line
    timer_end=$(date +%s)
    elapsed=$((timer_end - timer_start))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    runtime_line=$(printf "[RUNTIME] %02d:%02d:%02d" "$hours" "$minutes" "$seconds")
    if [[ -n "${log_file:-}" ]]; then
        echo "$runtime_line" >> "$log_file"
        if [[ "$status" -ne 0 ]]; then
            echo "[FAILED] $(date '+%F %T') exit=${status}" >> "$log_file"
        fi
    fi
    echo "$runtime_line"
    exit "$status"
}
trap finish EXIT

if [[ "$#" -ne 4 ]]; then
    echo "[ERROR] Usage: $0 <region_key> <final_tag> <sample_tag> <samples>" >&2
    exit 2
fi

region_key="$1"
final_tag="$2"
sample_tag="$3"
samples="$4"

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna}"
PLOT_DIR="${PROJECT_DIR}/Plot"
SCRIPTS_DIR="${PLOT_DIR}/scripts"
OUTPUT_DIR="${PLOT_DIR}/plots"
VARIABLES_DIR="${OUTPUT_DIR}/variables_dataVmc"
LOG_DIR="${OUTPUT_DIR}/logs_split"
SIDEBAND_REWEIGHT_JSON="${SIDEBAND_REWEIGHT_JSON:-${PROJECT_DIR}/HZaMVA/reweights/sideband_run3_iterative.json}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MAX_EVENTS="${MAX_EVENTS:-}"
SKIP_SYSTEMATICS="${SKIP_SYSTEMATICS:-0}"
DATA_VMC_BACKEND="${DATA_VMC_BACKEND:-auto}"
SETUP_CONDA_ENV="${SETUP_CONDA_ENV:-auto}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
LOCALIZE_CONDA_ENV="${LOCALIZE_CONDA_ENV:-auto}"
LOCAL_CONDA_DIR="${LOCAL_CONDA_DIR:-${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/higgs-alp-ana-conda}"
# conda env tarball：2026-09-17 從 AFS 搬到 /eos/project，改成 job 內 xrdcp 取回。
#
# 為什麼不放 AFS：1.1 GB 的單一大檔佔 AFS 配額（那裡只有 100 GB 且吃緊），
# 而 EOS 正好擅長單一大檔。⚠️ 反過來**不要**把解壓後的 env（十萬個小檔）放 EOS，
# tar 解壓到 EOS FUSE 會報 `Cannot close: Bad address` —— 規矩是「只放 tarball、
# 在節點本地解壓」。
#
# 為什麼是 /eos/project 而不是 /eos/cms/.../phys_susy：實測 2026-09-16，
# phys_susy 的 zh group 是 178.72/180.00 TB（99.29%、exceeded），
# /eos/project/h/htozg-dy-privatemc 是 11.82/20.00 TB（59.12%、ok）。
#
# 為什麼是 xrdcp 而不是直接 tar 讀 EOS 路徑：走 FUSE 讀在並行下比 AFS 還脆弱；
# xrdcp 先抓到節點本地 scratch 再解壓是 HZgamma B-mode 驗證過的作法
# （A-mode 直接讀 EOS 上的 env 實測失敗率 30.2%，B-mode 0%）。
CONDA_TARBALL_URL="${CONDA_TARBALL_URL:-root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/higgs-alp-ana.tar.gz}"
# 仍保留本地路徑的用法：設了 CONDA_TARBALL 就走它，不 xrdcp（離線／除錯用）。
CONDA_TARBALL="${CONDA_TARBALL:-}"

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
    source "$ANACONDA_SETUP"
    conda activate "$CONDA_ENV_NAME"
    set -u

    if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/python3" ]]; then
        echo "[ERROR] Failed to activate conda env: $CONDA_ENV_NAME" >&2
        exit 2
    fi

    PYTHON_BIN="${CONDA_PREFIX}/bin/python3"
    echo "[ENV] active CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
    echo "[ENV] active python=${PYTHON_BIN}"
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
        # 明確指定了本地路徑 -> 照舊直接讀（離線／除錯）
        if [[ ! -s "$CONDA_TARBALL" ]]; then
            echo "[ERROR] CONDA_TARBALL set but missing: $CONDA_TARBALL" >&2
            exit 2
        fi
        _pack="$CONDA_TARBALL"
        echo "[ENV] Extract conda tarball from $_pack to $LOCAL_CONDA_DIR"
    else
        # 預設：從 EOS xrdcp 到節點本地再解壓
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
    # 只刪自己抓下來的那份，別刪使用者指定的本地檔
    [[ -z "$CONDA_TARBALL" ]] && rm -f "$_pack"

    # 2026-09-23: export PATH BEFORE conda-unpack. conda-pack gives conda-unpack a
    # `#!/usr/bin/env python` shebang, so it needs a python on PATH. This only worked
    # here because SETUP_CONDA_ENV=auto had already activated the EOS env first; once
    # that activation is removed (see dataVmc.submit) there is no python and the job
    # dies with "/usr/bin/env: 'python': No such file or directory", return value 127.
    # The identical chain bit the p2root scoring jobs on 2026-09-22.
    export CONDA_PREFIX="$LOCAL_CONDA_DIR"
    export PATH="$LOCAL_CONDA_DIR/bin:$PATH"
    export LD_LIBRARY_PATH="$LOCAL_CONDA_DIR/lib:${LD_LIBRARY_PATH:-}"
    PYTHON_BIN="$LOCAL_CONDA_DIR/bin/python3"

    if [[ -x "$LOCAL_CONDA_DIR/bin/conda-unpack" ]]; then
        echo "[ENV] Run conda-unpack"
        "$LOCAL_CONDA_DIR/bin/conda-unpack"
    fi
}

mkdir -p "$LOG_DIR" "$VARIABLES_DIR"
export PYTHONPATH="${PYTHONPATH:-}:${PLOT_DIR}/lib:${PROJECT_DIR}/HZaMVA/scripts"
export PYTHONUNBUFFERED=1

partial_tag="${final_tag}_part_${sample_tag}"
log_file="${LOG_DIR}/${final_tag}_${region_key}_${sample_tag}.log"

cd "$PLOT_DIR"

{
    echo "[START] $(date '+%F %T')"
    echo "[HOST] $(hostname)"
    echo "[PWD] $(pwd)"
    echo "[JOB] tag=${final_tag} region=${region_key} sample_tag=${sample_tag} samples=${samples}"
    echo "[ENV] initial PYTHON_BIN=${PYTHON_BIN}"
    echo "[ENV] initial CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
    echo "[ENV] SETUP_CONDA_ENV=${SETUP_CONDA_ENV}"
    echo "[ENV] CONDA_ENV_NAME=${CONDA_ENV_NAME}"
    echo "[ENV] ANACONDA_SETUP=${ANACONDA_SETUP}"
    echo "[ENV] LOCALIZE_CONDA_ENV=${LOCALIZE_CONDA_ENV}"
    echo "[ENV] CONDA_TARBALL=${CONDA_TARBALL}"
    echo "[ENV] SKIP_SYSTEMATICS=${SKIP_SYSTEMATICS}"
    echo "[ENV] DATA_VMC_BACKEND=${DATA_VMC_BACKEND}"
} > "$log_file"

activate_conda_env_if_needed >> "$log_file" 2>&1
localize_conda_env_if_needed >> "$log_file" 2>&1
"$PYTHON_BIN" -u -c 'import numpy; import xgboost; from root_compat import import_pyroot; ROOT = import_pyroot(); ROOT.gROOT.SetBatch(True); print("[ENV] python imports OK")' >> "$log_file" 2>&1

cmd=(
    "$PYTHON_BIN" -u "$SCRIPTS_DIR/1_prepare_dataVmc.py"
    -y run3
    -m
    --ln
    --histOnly
    --backend "$DATA_VMC_BACKEND"
    --optimizeBranches
    --samples "$samples"
    --outputTag "$partial_tag"
)

case "$SKIP_SYSTEMATICS" in
    1|true|TRUE|yes|YES) cmd+=(--skipSystematics) ;;
    0|false|FALSE|no|NO) ;;
    *)
        echo "[ERROR] SKIP_SYSTEMATICS must be 0/1, true/false, or yes/no; got '$SKIP_SYSTEMATICS'" >&2
        exit 2
        ;;
esac

case "$region_key" in
    SR)  cmd+=(--region 1) ;;
    CR)  cmd+=(--region 2) ;;
    mva) cmd+=(-b) ;;
    *)
        echo "[ERROR] Unknown region key: $region_key" >&2
        exit 2
        ;;
esac

if [[ "$final_tag" == "sideband_rwgt" ]]; then
    cmd+=(--useSidebandReweight --sidebandReweightJson "$SIDEBAND_REWEIGHT_JSON" --noSidebandReweightUnc)
fi

if [[ -n "$MAX_EVENTS" ]]; then
    cmd+=(--maxEvents "$MAX_EVENTS")
fi

{
    printf "[CMD]"
    printf " %q" "${cmd[@]}"
    echo
} >> "$log_file"

"${cmd[@]}" >> "$log_file" 2>&1
echo "[DONE] $(date '+%F %T')" >> "$log_file"
