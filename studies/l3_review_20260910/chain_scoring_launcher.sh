#!/usr/bin/env bash
# Waits for the low/high-mass training to finish and PASS, then runs the scoring chain.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while [ "$(systemctl --user is-active hza-trainlh 2>/dev/null)" = "active" ]; do sleep 60; done
if ! grep -q "VERDICT: PASSED" $D/logs_fsrfix/train_lowhigh.log 2>/dev/null; then
  echo "[$(date '+%F %T')] low/high training did not pass -- scoring not started" \
    > $D/logs_fsrfix/scoring.log
  exit 1
fi
bash $D/chain_scoring.sh > $D/logs_fsrfix/scoring.log 2>&1
echo "RC=$?" >> $D/logs_fsrfix/scoring.log
