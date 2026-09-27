#!/usr/bin/env bash
#
# Re-merge the fpo=1 production after backfilling.
#
# WHY --force
#   The driver already wrote merged_nominal.parquet during the main run, before
#   the 81 backfilled jobs landed. Without --force the merger skips existing
#   merged files, which is exactly the merge-before-rechunk race: merge holds
#   only the old chunks, and the three-stage reconciliation would still pass
#   while being wrong. (CLAUDE.md: any backfill -> re-merge that sample.)
#
# WHY NOT 4_merge_parquet.sh
#   Its paths are hardcoded to parquet_DNA (the May 2026 production), the same
#   stale-path problem the p2root scripts have. Calling the underlying python
#   directly avoids editing a shared script for a one-off directory.
#
set -uo pipefail
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
N=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
L=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix
for st in Bkg_MC Data; do
    echo "=== merging $st at $(date '+%F %T') ==="
    ( cd "$R" && PYTHONPATH="$R" "$PY" scripts/A_merge_parquet_outputs.py "$N/$st" --force )
    echo "=== $st merger returned $? ==="
done
