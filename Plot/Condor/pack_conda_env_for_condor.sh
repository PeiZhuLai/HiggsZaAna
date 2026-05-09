#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANACONDA_SETUP="${ANACONDA_SETUP:-/eos/home-p/pelai/App/Anaconda/Anaconda/env_Anaconda.sh}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-higgs-alp-ana}"
ENV_CACHE_DIR="${ENV_CACHE_DIR:-${script_dir}/env_cache}"
CONDA_TARBALL="${CONDA_TARBALL:-${ENV_CACHE_DIR}/${CONDA_ENV_NAME}.tar.gz}"
FORCE_ENV_PACK="${FORCE_ENV_PACK:-0}"

if [[ -s "$CONDA_TARBALL" && "$FORCE_ENV_PACK" != "1" ]]; then
    echo "[ENV] Reuse existing tarball: $CONDA_TARBALL"
    echo "[ENV] Set FORCE_ENV_PACK=1 to rebuild it."
    exit 0
fi

if [[ ! -r "$ANACONDA_SETUP" ]]; then
    echo "[ERROR] Cannot read ANACONDA_SETUP: $ANACONDA_SETUP" >&2
    exit 2
fi

mkdir -p "$ENV_CACHE_DIR"

set +u
source "$ANACONDA_SETUP"
conda activate "$CONDA_ENV_NAME"
set -u

if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/python3" ]]; then
    echo "[ERROR] Failed to activate conda env: $CONDA_ENV_NAME" >&2
    exit 2
fi

tmp_tarball="${CONDA_TARBALL}.tmp"
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
