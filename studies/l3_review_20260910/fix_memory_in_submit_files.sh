#!/usr/bin/env bash
#
# Raise RequestMemory in the generated condor submit files, 3000 -> 12000 MB.
#
# WHY
#   The copy-then-read change (analysis.py stage_and_open) stages each NanoAOD to
#   node-local scratch before reading. Measured on the 20-job validation batch of
#   2026-09-15:
#       MemoryUsage  n=18  median 7325 MB  p90 9766 MB  max 9766 MB
#   against the old 3000 MB request -- which is why 18 of the 20 went straight to
#   a cgroup memory hold and only ran after release_memory_holds.sh bumped them.
#   Leaving 3000 in place for the full ~2000-job backfill would mean every job
#   burns one wasted run before the releaser rescues it.
#
#   12000 MB is p90 plus margin. The CPU-coercion trap was checked rather than
#   assumed: the bumped validation jobs came back
#       RequestMemory=8000 RequestCpus=1 CpusProvisioned=1 MemoryProvisioned=8000
#   so this pool does not silently widen a single-core request here. Re-verify
#   after submitting anyway -- a 2-core coercion doubles the fair-share cost and
#   makes slots much harder to match.
#
# WHY NOT --reconfigure_jobs
#   It does not rewrite submit files. managers.py:470 calls job.submit() without
#   reconfigure=True, so jobs.py:240 skips the writer when the files exist. The
#   log claims otherwise. Same reason the JobFlavour fix had to be done this way.
#
# USAGE
#   bash fix_memory_in_submit_files.sh             # dry run: count only
#   DRY_RUN=0 bash fix_memory_in_submit_files.sh   # rewrite
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
NEW="${NEW:-12000}"
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
OUT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/fix_memory.log

mapfile -t subs < <(find "$R/eos_logs/Bkg_MC" "$R/eos_logs/Data" -name '*_batch_submit_job*.txt' 2>/dev/null)
echo "submit files found: ${#subs[@]}"
if [ "$DRY_RUN" = "1" ]; then
    echo "current RequestMemory distribution:"
    grep -h 'RequestMemory' "${subs[@]}" 2>/dev/null | sort | uniq -c
    echo "would set RequestMemory = ${NEW} in all of them"
    exit 0
fi
# RAISE ONLY. 105 of these files are already at 20000 because
# release_memory_holds.sh climbed the ladder for jobs that genuinely needed it;
# a blanket sed would quietly demote them back to 12000 and they would hold
# again on their next run. Touch a file only when its current value is lower.
raised=0; kept=0
for f in "${subs[@]}"; do
    cur=$(sed -n 's/^RequestMemory = \([0-9]*\).*/\1/p' "$f" | head -1)
    if [ -n "$cur" ] && [ "$cur" -lt "$NEW" ] 2>/dev/null; then
        sed -i "s/^RequestMemory = .*/RequestMemory = ${NEW}/" "$f" && raised=$((raised+1))
    else
        kept=$((kept+1))
    fi
done
echo "$(date '+%F %T') raised ${raised} submit files to ${NEW} MB, left ${kept} at their existing (>=) value" | tee -a "$OUT"
echo "verification straight off disk:" | tee -a "$OUT"
grep -h 'RequestMemory' "${subs[@]}" 2>/dev/null | sort | uniq -c | tee -a "$OUT"
