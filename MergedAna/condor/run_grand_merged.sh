#!/bin/bash
# condor payload for the merged limit chain (1_grand_merged.sh).
# Stage switches (RUN_*) and MASSES arrive through the .sub `environment`.
# The chain writes its products to AFS/EOS itself. Its log goes to AFS one line
# per open/append/close: AFS publishes writes to other clients only on close,
# so a long-lived `tee` stays invisible from lxplus until the job ends.
set -o pipefail
CHAIN_LOG="${CHAIN_LOG:?CHAIN_LOG not set}"
logl() { while IFS= read -r l; do printf '%s\n' "$l" >> "$CHAIN_LOG"; done; }
echo "host=$(hostname) start=$(date '+%F %T')" > "$CHAIN_LOG"
env | grep -E '^(RUN_|MASSES=|MERGED_)' | sort | logl
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna || exit 2
bash 1_grand_merged.sh 2>&1 | logl
rc=${PIPESTATUS[0]}
echo "chain exit code: ${rc}  end=$(date '+%F %T')" >> "$CHAIN_LOG"
exit "$rc"
