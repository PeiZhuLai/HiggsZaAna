#!/bin/bash
# Formal/production NN training: per-mass DisCo (disco_masses=1,2,3) at lambda=100, the chosen
# best setting. epochs=30 -- the VALIDATED config (job 11642604 -> AUC 0.9806, R(mA=1)=1.50).
# 40 epochs drifts R(mA=1) up to ~2.0 (weak DisCo => run-to-run variance), so fix epochs=30 and
# scan seeds; a picker selects the lowest-R(mA=1) model afterwards.
# Arg: SEED. Writes per-seed model + loss/overtraining PDFs.
#
# ROBUST conda activation: the hza_NN env lives on EOS; activation + first import on a fresh
# worker is flaky (intermittent OSError errno5 / "No module named torch" when EOS isn't warm).
# Retry activate + import a few times, and FAIL LOUDLY (exit 1) instead of a silent "done".
set -o pipefail
SEED="${1:-1}"
echo "==== FINAL NN train (per-mass lam100, epochs=30, seed=${SEED}) on $(hostname) at $(date '+%F %T') ===="

source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh

ok=0
for attempt in 1 2 3 4 5 6; do
    conda activate hza_NN 2>/dev/null
    if python3 -c "import numpy, uproot, torch, torch.nn" 2>/dev/null; then
        ok=1; break
    fi
    echo "[retry ${attempt}] hza_NN not importable yet (EOS not warm); sleeping 20s..."
    conda deactivate 2>/dev/null; sleep 20
done
[ "$ok" = 1 ] || { echo "[FATAL] could not import numpy/uproot/torch from hza_NN after retries"; exit 1; }
echo "[env] $CONDA_DEFAULT_ENV  python=$(python3 --version 2>&1)"

echo "==== GPU ===="
nvidia-smi || echo "[WARN] nvidia-smi unavailable -- not on a GPU node?"
python3 -c "import torch; print('[torch]', torch.__version__, 'cuda available:', torch.cuda.is_available(), '|', (torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no GPU'))"

cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts
OUT="model_Za_NN_run3_permass_lam100_s${SEED}.pt"

echo "==== train (per-mass DisCo lam=100, disco_masses=1,2,3, full sample, 30 epochs, seed=${SEED}) ===="
python3 train_nn.py --lambda-disco 100 --epochs 30 --seed "${SEED}" \
    --n-sig-per 60000 --n-bkg 740000 --disco-masses 1,2,3 \
    --out "${OUT}" --plot-dir plots_nn
rc=$?

[ $rc -eq 0 ] && [ -f "${OUT}" ] || { echo "[FATAL] training failed (rc=$rc) or model ${OUT} missing"; exit 1; }
echo "==== done seed=${SEED} (model ${OUT}) at $(date '+%F %T') ===="
