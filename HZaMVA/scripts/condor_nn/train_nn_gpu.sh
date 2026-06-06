#!/bin/bash
# Executable for the GPU condor job. Sets up hza_NN, logs the GPU, runs NN training.
# NOT set -u (conda activate references unbound vars).
set -o pipefail
LAMBDA="${1:-0}"          # DisCo strength, passed from the submit file (queue LAMBDA in 0 50)
echo "==== NN train lambda=${LAMBDA} on $(hostname) at $(date '+%F %T') ===="

source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
conda activate hza_NN
echo "[env] $CONDA_DEFAULT_ENV  python=$(python3 --version 2>&1)"

echo "==== GPU ===="
nvidia-smi || echo "[WARN] nvidia-smi unavailable -- not on a GPU node?"
python3 -c "import torch; print('[torch]', torch.__version__, 'cuda available:', torch.cuda.is_available(), '|', (torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no GPU'))"

cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts

echo "==== train (lambda-disco=${LAMBDA}, full sample, 30 epochs) ===="
# parametrized MLP + DisCo + sideband reweight, GPU. ~60k signal/mass + full background.
python3 train_nn.py --lambda-disco "${LAMBDA}" --epochs 40 \
    --n-sig-per 60000 --n-bkg 740000 --disco-masses 1,2,3 \
    --out model_Za_NN_run3_permassXL_lam${LAMBDA}.pt

echo "==== done lambda=${LAMBDA} at $(date '+%F %T') ===="
