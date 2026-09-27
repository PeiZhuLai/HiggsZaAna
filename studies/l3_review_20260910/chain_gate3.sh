#!/usr/bin/env bash
# Waits for hza-p2root to finish, then runs GATE 3 automatically.
# Runs as its own systemd unit so it survives logout (el9 KillUserProcesses=yes).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while [ "$(systemctl --user is-active hza-p2root 2>/dev/null)" = "active" ]; do sleep 120; done
echo "[$(date '+%F %T')] hza-p2root finished, running GATE 3" >> $D/logs_fsrfix/chain.log
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
python $D/gate3_reconcile.py --since 1789906080 >> $D/logs_fsrfix/chain.log 2>&1
echo "[$(date '+%F %T')] GATE 3 rc=$?" >> $D/logs_fsrfix/chain.log
