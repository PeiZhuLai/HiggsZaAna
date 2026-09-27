#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ENV_CACHE_DIR="${ENV_CACHE_DIR:-${script_dir}/env_cache}"
CONDA_TARBALL="${CONDA_TARBALL:-${ENV_CACHE_DIR}/${CONDA_ENV_NAME}.tar.gz}"
FORCE_ENV_PACK="${FORCE_ENV_PACK:-0}"

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

# 快取判斷：tarball 現在住在 EOS，所以要問 EOS 而不是本地路徑。
# 🔴 這裡一定要跟著改：上傳成功後本地檔會被刪掉，如果還用 `-s "$CONDA_TARBALL"`
# 判斷，本地永遠不存在 -> 每次提交都重打包一次 1.1 GB 的 env，快取等於失效
# （而且不會有任何錯誤訊息，只是每次都慢）。
CONDA_TARBALL_URL="${CONDA_TARBALL_URL:-root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/${CONDA_ENV_NAME}.tar.gz}"
_eos_has_tarball() { xrdfs "${CONDA_TARBALL_URL%%//eos/*}" stat "/eos/${CONDA_TARBALL_URL#*//eos/}" >/dev/null 2>&1; }

if [[ "$FORCE_ENV_PACK" != "1" ]] && { [[ -s "$CONDA_TARBALL" ]] || _eos_has_tarball; }; then
    if [[ -s "$CONDA_TARBALL" ]]; then
        echo "[ENV] Reuse existing tarball: $CONDA_TARBALL"
    else
        echo "[ENV] Reuse existing tarball on EOS: $CONDA_TARBALL_URL"
    fi
    echo "[ENV] Set FORCE_ENV_PACK=1 to rebuild it."
    exit 0
fi

if [[ ! -r "$ANACONDA_SETUP" ]]; then
    echo "[ERROR] Cannot read ANACONDA_SETUP: $ANACONDA_SETUP" >&2
    exit 2
fi

mkdir -p "$ENV_CACHE_DIR"

set +u
sanitize_python_env_for_conda
source "$ANACONDA_SETUP"
conda activate "$CONDA_ENV_NAME"
set -u

if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/python3" ]]; then
    echo "[ERROR] Failed to activate conda env: $CONDA_ENV_NAME" >&2
    exit 2
fi

case "$CONDA_TARBALL" in
    *.tar.gz) tmp_tarball="${CONDA_TARBALL%.tar.gz}.tmp.tar.gz" ;;
    *.tgz)    tmp_tarball="${CONDA_TARBALL%.tgz}.tmp.tgz" ;;
    *)        tmp_tarball="${CONDA_TARBALL}.tmp.tar.gz" ;;
esac
rm -f "$tmp_tarball"

echo "[ENV] Pack conda env: $CONDA_PREFIX"
echo "[ENV] Output tarball: $CONDA_TARBALL"

if command -v conda-pack >/dev/null 2>&1; then
    conda-pack -p "$CONDA_PREFIX" -o "$tmp_tarball" --force
elif "$CONDA_PREFIX/bin/python3" -c 'import conda_pack' >/dev/null 2>&1; then
    "$CONDA_PREFIX/bin/python3" -m conda_pack -p "$CONDA_PREFIX" -o "$tmp_tarball" --force
else
    echo "[WARN] conda-pack is not available; using plain tar fallback." >&2
    echo "[WARN] If jobs fail after extraction, install conda-pack and rebuild this tarball." >&2
    tar -czf "$tmp_tarball" -C "$CONDA_PREFIX" .
fi

mv "$tmp_tarball" "$CONDA_TARBALL"
ls -lh "$CONDA_TARBALL"

# ── 上傳到 EOS，並把 AFS 上的本地檔刪掉 ──────────────────────────────────────
# 2026-09-17：job 端（run_dataVmc_condor_job.sh）已改成從 EOS xrdcp 取回，
# 所以打包完要送上去，否則下次重打包又會在 AFS 留下 1.1 GB。
# 為什麼是 /eos/project 而不是 /eos/cms/.../phys_susy：後者的 zh group 實測
# 178.72/180.00 TB（99.29%、exceeded），前者 11.82/20.00 TB（59.12%、ok）。
# 設 SKIP_EOS_UPLOAD=1 可跳過（純本地除錯）。
CONDA_TARBALL_URL="${CONDA_TARBALL_URL:-root://eosproject-h.cern.ch//eos/project/h/htozg-dy-privatemc/pelai/App/${CONDA_ENV_NAME}.tar.gz}"
if [[ "${SKIP_EOS_UPLOAD:-0}" != "1" ]]; then
    echo "[ENV] Upload to $CONDA_TARBALL_URL"
    if xrdcp -f "$CONDA_TARBALL" "$CONDA_TARBALL_URL"; then
        # 驗證再刪：只比大小不夠，壞掉的複本大小照樣相同，所以抓回來 gzip -t。
        _v="${TMPDIR:-/tmp}/.verify_$$.tar.gz"
        if xrdcp -f -s "$CONDA_TARBALL_URL" "$_v" \
           && [[ "$(stat -c %s "$_v")" == "$(stat -c %s "$CONDA_TARBALL")" ]] \
           && gzip -t "$_v" 2>/dev/null; then
            rm -f "$_v"
            echo "[ENV] EOS copy verified (size + gzip -t); removing local $CONDA_TARBALL"
            rm -f "$CONDA_TARBALL"
        else
            rm -f "$_v"
            echo "[WARN] EOS copy failed verification -- keeping local tarball at $CONDA_TARBALL" >&2
        fi
    else
        echo "[WARN] xrdcp upload failed -- keeping local tarball at $CONDA_TARBALL" >&2
    fi
fi
