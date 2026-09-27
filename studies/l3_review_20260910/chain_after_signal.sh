#!/usr/bin/env bash
#
# Wait for the signal re-production to finish, reconcile it, and only then
# submit the bkg + data stages.
#
# WHY the gate is fail-closed:
#   run_analysis.py returns 0 even when every job failed -- it counts jobs that
#   were retired after 5 worker failures as "COMPLETED : 100.00 percent". The
#   first submission attempt today returned 0 with 1104/1104 jobs dead (exit
#   127) and zero output. So the driver's exit status proves nothing; the gate
#   is the artifact count plus the chunk==merged reconciliation.
#
# Gate (all must pass before 6762 more jobs are submitted):
#   1. the signal driver process is gone;
#   2. >= MIN_PARQUET nominal chunks exist;
#   3. reconcile_signal.py exits 0 (chunk==merged, merged not older than its
#      newest chunk, and <1% drift against the previous production).
#
# USAGE
#   bash chain_after_signal.sh            # dry run: waits, reconciles, does NOT submit
#   SUBMIT=1 bash chain_after_signal.sh   # submit bkg+data if the gate passes
#
set -uo pipefail

D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910"
STAG="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC"
PY="/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python"
SUBMIT="${SUBMIT:-0}"
MIN_PARQUET="${MIN_PARQUET:-1090}"     # of 1104; a handful of transient xrootd
                                        # losses is tolerable, a systematic hole is not
POLL="${POLL:-300}"

echo "[chain] waiting for the signal driver to finish (poll ${POLL}s)"
while ps -u "$(id -un)" -o cmd= 2>/dev/null | grep -q "[r]un_analysis.py --config metadata/za_signal_run3.json"; do
    sleep "$POLL"
done
echo "[chain] signal driver has exited at $(date '+%F %T')"

n=$(find "$STAG" -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l | tr -d ' ')
echo "[chain] nominal chunks on disk: ${n}/1104"
if [ "${n:-0}" -lt "$MIN_PARQUET" ]; then
    echo "[chain] GATE FAILED: only ${n} chunks (need >= ${MIN_PARQUET}). Not submitting."
    echo "[chain] Likely retired jobs -- inspect eos_logs/Sig_MC/*/*/*.err before rerunning."
    exit 1
fi

echo "[chain] running three-stage reconciliation"
"$PY" "$D/reconcile_signal.py" | tee "$D/logs_fsrfix/reconcile_signal.txt"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then
    echo "[chain] GATE FAILED: reconciliation exited ${rc}. Not submitting."
    exit 1
fi
echo "[chain] GATE PASSED"

# Direct resolution measurement, before the long downstream chain. The electron
# channel is the control: the FSR recovery is muon-only, so it must not move.
echo "[chain] measuring sigma_eff, old vs new"
"$PY" "$D/compare_sigma_eff.py" > "$D/logs_fsrfix/sigma_eff_stdout.txt" 2>&1
tail -6 "$D/logs_fsrfix/sigma_eff_stdout.txt" | sed 's/^/[chain] /'

if [ "$SUBMIT" != "1" ]; then
    echo "[chain] SUBMIT is not 1 -- stopping here without submitting."
    exit 0
fi

# bkg+data go through the wave submitter, not submit_fsr_reproduction.sh: the
# signal stage showed that ~1089 concurrent jobs make the shared conda env on
# /eos return I/O errors for 66 % of attempts, while ~211 concurrent was clean.
# The wave script raises fpo (6762 -> ~981 jobs, so 4-8x fewer conda-env reads)
# and runs one (stage, year) wave at a time, each sized near 250-300 jobs.
echo "[chain] handing over to the wave submitter at $(date '+%F %T')"
DRY_RUN=0 bash "$D/submit_bkg_data_waves.sh" >> "$D/logs_fsrfix/waves.log" 2>&1
echo "[chain] wave submitter returned $? (NOT a success indicator -- check artifacts)"
echo "[chain] bkg+data finished at $(date '+%F %T')"
