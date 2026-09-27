#!/usr/bin/env bash
#
# Wait for the fpo=1 production, reconcile it, and only then convert to ROOT.
#
# It deliberately STOPS after p2root. Everything downstream of that (BDT apply,
# fTest, signalFit, plotEffSigma) changes analysis results, and this week has
# produced enough evidence that an unattended chain will happily carry a broken
# input all the way to a plot: run_analysis.py returned 0 with every job dead,
# --reconfigure_jobs logged a rewrite it never did, and a driver reported
# "100.00 percent completed" while 2028 chunks were missing. So the gates are
# fail-closed and the last step is a report, not a fit.
#
# GATES
#   1. registry says the production finished (watchdog_fpo1 writes DONE/STOPPED)
#   2. reconcile_fpo1_all.sh exits 0
#        rc=1 production incomplete (job dirs vs chunks)
#        rc=2 pipeline loss (chunk_rows != merged_rows, or merge-before-rechunk)
#   3. after p2root, root events == merged rows for every sample
#
# USAGE
#   bash chain_after_fpo1.sh            # dry run: shows what it would do
#   RUN=1 nohup bash chain_after_fpo1.sh >/dev/null 2>&1 &
#
set -uo pipefail
RUN="${RUN:-0}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
REG=$D/logs_fsrfix/PENDING_RESULTS.txt
LOG=$D/logs_fsrfix/chain_fpo1.log

say(){ echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

if [ "$RUN" != "1" ]; then
    echo "would: wait for 'DONE  fpo1' in $REG"
    echo "  then: bash $D/reconcile_fpo1_all.sh   (gate, must exit 0)"
    echo "  then: bash $REPO/1_run_P2Root.sh && bash $REPO/2_prepare_rootfile.sh"
    echo "  then: verify root events == merged rows, write result to $REG"
    echo "  and STOP -- fitting is left to a human decision"
    exit 0
fi

say "waiting for the production to finish"
while ! grep -qE 'DONE  fpo1|STOPPED  fpo1' "$REG" 2>/dev/null; do sleep 300; done
if grep -q 'STOPPED  fpo1' "$REG"; then
    say "GATE 1 FAILED: production stopped short -- see $REG"; exit 1
fi
say "GATE 1 passed"

say "running reconciliation"
bash "$D/reconcile_fpo1_all.sh" >/dev/null 2>&1
rc=$?
if [ "$rc" -ne 0 ]; then
    say "GATE 2 FAILED: reconcile rc=$rc (1=incomplete, 2=pipeline loss). Not converting."
    echo "$(date '+%F %H:%M') BLOCKED fpo1-chain  reconcile rc=$rc" >> "$REG"
    exit 2
fi
say "GATE 2 passed"

say "p2root"
( cd "$REPO" && bash 1_run_P2Root.sh && bash 2_prepare_rootfile.sh ) >> "$LOG" 2>&1
say "p2root returned $? (not a success indicator -- the artifact check follows)"
echo "$(date '+%F %H:%M') DONE  fpo1-chain  p2root finished, root-vs-merged check pending" >> "$REG"
