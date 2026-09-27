#!/usr/bin/env bash
#
# Rewrite +JobFlavour in the already-generated condor submit files, longlunch -> workday.
#
# WHY NOT --reconfigure_jobs
#   It does not do this. run_analysis.py --reconfigure_jobs reaches
#   AnalysisManager.modify_jobs (analysis.py:427), which logs
#   "Forcing reconfiguration of jobs (rewriting python configs, executables and
#   condor_submit files)" and calls submit_jobs(dry_run=True) with
#   customized_jobs=False. But customize_jobs() only recomputes PATHS, and the
#   actual writer is reached through
#       managers.py:470  job.submit(dry_run = True)
#   which never passes reconfigure=True. jobs.py:240 therefore takes its other
#   branch -- the files already exist, so nothing is rewritten. Verified: after a
#   full reconfigure run the submit files were still stamped Sep 12 20:37 and
#   condor_q reported MaxRuntime=7200 / longlunch for all 6071 jobs.
#   The log line says it rewrote them. It did not. (Verify the artifact.)
#
# WHY IT MATTERS
#   Walltime is fixed at SUBMIT time; condor_qedit on a live job is ignored. On
#   the night of 2026-09-13, of 6000 HZa jobs 5081 were removed for wall time and
#   only 919 finished -- and those that finished needed a median of 5600 s with a
#   p90 of 10495 s, against a 7200 s wall. The job-time distribution straddles the
#   limit, so the 8 h flavour should convert most of the losses into completions.
#   Jobs removed for wall time return NO stdout/stderr, so this failure mode is
#   invisible except in condor history.
#
# USAGE
#   bash fix_walltime_in_submit_files.sh             # dry run: count only
#   DRY_RUN=0 bash fix_walltime_in_submit_files.sh   # rewrite
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
OUT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/fix_walltime.log

# Scoped to the two stage trees. eos_logs as a whole holds hundreds of thousands
# of files from old productions and cannot be walked in reasonable time.
mapfile -t subs < <(find "$R/eos_logs/Bkg_MC" "$R/eos_logs/Data" -name '*_batch_submit_job*.txt' 2>/dev/null)
echo "submit files found: ${#subs[@]}"
n_long=0
for f in "${subs[@]}"; do grep -q 'longlunch' "$f" 2>/dev/null && n_long=$((n_long+1)); done
echo "carrying longlunch : ${n_long}"

if [ "$DRY_RUN" = "1" ]; then
    echo "would rewrite ${n_long} files: +JobFlavour = \"longlunch\" -> \"workday\""
    exit 0
fi
changed=0
for f in "${subs[@]}"; do
    if grep -q 'longlunch' "$f" 2>/dev/null; then
        sed -i 's/+JobFlavour = "longlunch"/+JobFlavour = "workday"/' "$f" && changed=$((changed+1))
    fi
done
echo "$(date '+%F %T') rewrote ${changed} submit files longlunch -> workday" | tee -a "$OUT"
# Prove it on disk rather than trusting the loop.
left=$(grep -l 'longlunch' "${subs[@]}" 2>/dev/null | wc -l)
echo "$(date '+%F %T') submit files still carrying longlunch: ${left} (must be 0)" | tee -a "$OUT"
