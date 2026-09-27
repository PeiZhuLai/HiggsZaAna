#!/usr/bin/env bash
# Runs GATE 4b (branch closure) once the hadd stage finishes. Separate unit so it does
# not race the GATE 4 call that lives inside run_prepare_fsrfix.sh.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while [ "$(systemctl --user is-active hza-prepare 2>/dev/null)" = "active" ]; do sleep 60; done
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
python $D/gate4b_branch_closure.py > $D/logs_fsrfix/gate4b.log 2>&1
echo "exit=$?" >> $D/logs_fsrfix/gate4b.log
