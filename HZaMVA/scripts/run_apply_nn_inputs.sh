#!/bin/bash
# Step (b): apply the trained NN to the aggregated bdt_inputs files (the only files that carry
# the full _oHm feature set + inclusive/train/validation/test trees) and write slim, friend-style
# NN-scored copies for diagnostics (ROC / sculpting / overtraining vs the BDT).
# NOT the production analysis files -- that is step (a) (integrate into Parque2Root_BDT.py).
#
#   ./run_apply_nn_inputs.sh [MODEL.pt]
set -o pipefail
MODEL="${1:-/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_NN_run3_permass_lam100.pt}"
INBASE="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
OUTBASE="/eos/home-p/pelai/HZa/root_P2Root/run3_nn_scored_inputs_nominal"
SCRIPTS="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts"

source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
conda activate hza_NN
cd "$SCRIPTS"
echo "[cfg] model=$MODEL"
[ -f "$MODEL" ] || { echo "[ERR] model not found: $MODEL"; exit 1; }

mkdir -p "$OUTBASE"
n=0; fail=0
for f in "$INBASE"/*/run3.root; do          # only dirs that actually carry features have run3.root
    d=$(basename "$(dirname "$f")")
    out="$OUTBASE/$d/run3.root"
    mkdir -p "$OUTBASE/$d"
    echo "==== [$((n+1))] $d ===="
    if python3 apply_nn.py --model "$MODEL" --input "$f" --output "$out" --slim --device cpu; then
        n=$((n+1))
    else
        echo "[ERR] failed: $d"; fail=$((fail+1))
    fi
done
echo "==== done: $n ok, $fail failed -> $OUTBASE ===="
