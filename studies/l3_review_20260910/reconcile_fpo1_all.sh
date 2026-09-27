#!/usr/bin/env bash
#
# Post-production reconciliation for the fpo=1 re-production.
#
# Two DIFFERENT questions, both required (CLAUDE.md):
#   A. completeness  -- did production make everything it was supposed to?
#      job dirs vs chunks on disk, per sample. chunk==merged says nothing about
#      this: a stage that lost 30 % of its jobs reconciles perfectly.
#   B. pipeline loss -- chunk_rows == merged_rows, plus the merge-before-rechunk
#      guard (merged mtime >= newest chunk mtime). reconcile_signal.py does this
#      and is already generic over --base/--ref.
# The root stage (== root_events) is checked after p2root, not here.
#
# Reference for the drift column is the PREVIOUS production, not the fpo=4
# attempt: that one is incomplete by construction (5059 of 6762 chunks) and
# comparing against it would report a fake deficit everywhere.
#
# USAGE
#   bash reconcile_fpo1_all.sh          # read-only, prints and writes a report
#
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
NEW=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1
OLD=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp
OLDDATA=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
REP=$D/logs_fsrfix/reconcile_fpo1.txt
REG=$D/logs_fsrfix/PENDING_RESULTS.txt

{
echo "=== fpo=1 reconciliation  $(date '+%F %T') ==="
rc=0
for st in Bkg_MC Data; do
    echo
    echo "############ $st : A. completeness (job dirs vs chunks) ############"
    $PY - "$NEW/$st" <<'PYEOF'
import os, sys, glob
base = sys.argv[1]
tj = tc = 0
for smp in sorted(os.listdir(base)):
    d = os.path.join(base, smp)
    if not os.path.isdir(d): continue
    # Completion marker is the SUMMARY json, not the parquet. At fpo=1 a job
    # covers one input file and "zero events selected" is routine -- HiggsDNA
    # writes the summary and no parquet, and the job is successful. Counting
    # parquets reported 2286 healthy jobs as missing on 2026-09-16 and triggered
    # 200 needless resubmissions. (ref_progress_count_needs_completion_marker)
    jobs = sorted(x for x in os.listdir(d) if x.startswith("job_"))
    notrun = [j for j in jobs
              if not glob.glob(os.path.join(d, j, "*_summary_%s.json" % j.replace("job_", "job")))]
    chunks = len([f for f in glob.glob(os.path.join(d, "**", "*_nominal.parquet"), recursive=True)
                  if os.path.basename(f) != "merged_nominal.parquet"])
    tj += len(jobs); tc += len(jobs) - len(notrun)
    flag = "" if not notrun else "   NOT RUN %d" % len(notrun)
    print("  %-34s jobs=%6d ran=%6d parquet=%6d%s"
          % (smp, len(jobs), len(jobs) - len(notrun), chunks, flag))
print("  TOTAL jobs=%d ran=%d not-run=%d" % (tj, tc, tj - tc))
sys.exit(1 if tj - tc else 0)
PYEOF
    [ $? -ne 0 ] && { echo "  [WARN] $st incomplete -- see MISSING lines above"; rc=1; }

    echo
    echo "############ $st : B. chunk==merged + mtime guard ############"
    # The drift reference differs per stage: parquet_DNA_tmp has Bkg_MC and
    # Sig_MC but NO Data (checked 2026-09-18), so Data must be compared against
    # parquet_DNA, which holds the only Data merged parquet of the pre-FSR-fix
    # production. Pointing both at parquet_DNA_tmp would silently report every
    # Data sample as having no reference.
    ref="$OLD/$st"; [ "$st" = "Data" ] && ref="$OLDDATA/$st"
    $PY "$D/reconcile_signal.py" --base "$NEW/$st" --ref "$ref" --max-drift 2.0
    [ $? -ne 0 ] && { echo "  [FAIL] $st pipeline reconciliation failed"; rc=2; }
done
echo
echo "=== overall rc=$rc (0 clean, 1 incomplete production, 2 pipeline loss) ==="
exit $rc
} 2>&1 | tee "$REP"
rc=${PIPESTATUS[0]}
echo "$(date '+%F %H:%M') RECONCILE fpo1 rc=$rc -> $REP" >> "$REG"
exit $rc
