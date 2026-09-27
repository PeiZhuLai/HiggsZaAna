#!/usr/bin/env bash
# Waits for the DY 2022 backfill, then re-runs GATE 3 over the full 1207 expected set.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while [ "$(systemctl --user is-active hza-dy2022 2>/dev/null)" = "active" ]; do sleep 60; done
echo "[$(date '+%F %T')] backfill done, re-running GATE 3" >> $D/logs_fsrfix/chain.log
set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
python $D/gate3_reconcile.py --since $(date -d "2026-09-20 14:08:00" +%s) >> $D/logs_fsrfix/chain.log 2>&1
echo "[$(date '+%F %T')] GATE 3 rerun rc=$?" >> $D/logs_fsrfix/chain.log
