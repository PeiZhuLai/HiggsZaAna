#!/bin/bash

# Set proxy
export X509_USER_PROXY=$PWD/GRID_PROXY

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/etc/profile.d/conda.sh" ]; then
        . "/usr/etc/profile.d/conda.sh"
    else
        export PATH="/usr/bin:$PATH"
    fi
fi
unset __conda_setup

# Keep ROOT/cling isolated from user-site Python packages and stray include paths.
export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset CPATH
unset CPLUS_INCLUDE_PATH
unset C_INCLUDE_PATH
export PYTHONPATH="HIGGSDNA_BASE"

# --- B-mode: if HZA_BMODE_PACK is set, fetch the conda-pack tarball to node-local
# scratch, extract, and activate it (avoids many workers reading the eos-fuse conda
# env concurrently). Kept last so the activated env's python wins PATH. Opt-in only:
# when unset, the classic eos-path python (jobs.py rewrite) is used unchanged.
if [ -n "$HZA_BMODE_PACK" ]; then
    ENVROOT="${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/hza_ana_env"
    mkdir -p "$ENVROOT"
    for _t in 1 2 3; do
        xrdcp -f -s "$HZA_BMODE_PACK" "$ENVROOT/pack.tar.gz" && break
        sleep 15
    done
    [ -s "$ENVROOT/pack.tar.gz" ] || { echo "ENV_FETCH_FAIL"; exit 110; }
    tar xzf "$ENVROOT/pack.tar.gz" -C "$ENVROOT" || { echo "ENV_EXTRACT_FAIL"; exit 111; }
    rm -f "$ENVROOT/pack.tar.gz"
    source "$ENVROOT/bin/activate"
    "$ENVROOT/bin/conda-unpack" 2>/dev/null || true
    export PYTHONPATH="HIGGSDNA_BASE"
fi

python -s PYTHON_FILE
