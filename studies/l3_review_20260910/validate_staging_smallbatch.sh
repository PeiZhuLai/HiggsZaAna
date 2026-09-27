#!/usr/bin/env bash
#
# Small-batch validation of the copy-then-read change before any mass resubmit.
#
# WHY A SMALL BATCH
#   The change is in the shared read path (analysis.py stage_and_open), so any
#   job exercises it. Resubmitting all ~2000 missing jobs to find out would cost
#   another night and another round of fair-share debt if it regresses. 20 jobs
#   from the sample with the worst deficit answers the same question in under an
#   hour.
#
# WHY THIS SAMPLE
#   DYJetsToLL_2022postEE is missing 186 of 421 chunks and is the one that
#   failed in the read probe with exactly the production failure:
#     STREAM FAIL 122.3 s  "File did not vector_read properly: Operation expired"
#     XRDCP  ok     6.5 s  -> local read 1.6 s
#
# WHY IT SUBMITS THE EXISTING .sub FILES DIRECTLY
#   Going through run_ana_bkgmc.sh would not restrict the batch: with
#   CLEAN_ANALYSIS_STATE=0 the pickled sample_list wins over the CLI
#   ("kwarg 'sample_list' ... will be ignored in favor of the pickled value"),
#   so the driver would submit every missing job in the stage. The per-job .sub
#   files already exist, already carry workday, and PYTHONPATH points at the
#   repo, so they pick up the new code with no regeneration.
#
# JUDGEMENT (decided before looking)
#   Baseline for these jobs, measured 2026-09-15 over 250 recent exit-1s:
#     vector_read timeout 237 | file unavailable 12 | other 1
#   and over all recent 8 h jobs: exit0 785 / exit1 633 / WALL 386 = 56 % loss.
#   PASS  = at least 15 of 20 produce a chunk, and no exit-1 carries
#           "vector_read properly".
#   FAIL  = vector_read failures persist -> the staging change is not the fix,
#           stop and reconsider (IHEP, or fewer files per job).
#
# USAGE
#   bash validate_staging_smallbatch.sh             # dry run
#   DRY_RUN=0 bash validate_staging_smallbatch.sh   # submit 20 jobs
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
N="${N:-20}"
SAMPLE="${SAMPLE:-DYJetsToLL_2022postEE}"
STAGE="${STAGE:-Bkg_MC}"
S=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/validate_staging.log

missing=()
for jd in "$S/$STAGE/$SAMPLE"/job_*; do
    [ -d "$jd" ] || continue
    n=$(basename "$jd"); n=${n#job_}
    [ -f "$jd/output_job_${n}_nominal.parquet" ] && continue
    sub="$R/eos_logs/$STAGE/$SAMPLE/job_$n/${SAMPLE}_batch_submit_job${n}.txt"
    [ -f "$sub" ] && missing+=("$sub")
    [ "${#missing[@]}" -ge "$N" ] && break
done
echo "sample=$SAMPLE  jobs without a chunk, selected: ${#missing[@]}"
[ "${#missing[@]}" -eq 0 ] && { echo "nothing to validate"; exit 1; }
echo "first: ${missing[0]}"
grep -h -E 'JobFlavour|RequestMemory' "${missing[0]}"

if [ "$DRY_RUN" = "1" ]; then
    echo "would condor_submit ${#missing[@]} existing .sub files (batch name left as-is)"
    exit 0
fi
# submit_template.txt is getenv = True, so the job inherits the SUBMITTING
# shell's environment. Submitting these .sub files from a bare shell gave all 20
# jobs exit 127 in 3 s: the wrapper leaves "python" literal (exe_template
# activates the env node-locally) and without HZA_BMODE_PACK it never activates
# anything, so there is no python on PATH. The driver normally exports these;
# doing it by hand means exporting them by hand too.
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$R"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"
export HZA_STAGE_INPUTS=1
ok=0
for sub in "${missing[@]}"; do
    condor_submit "$sub" >/dev/null 2>&1 && ok=$((ok+1))
done
echo "$(date '+%F %T') submitted ${ok}/${#missing[@]} validation jobs for ${SAMPLE}" | tee -a "$LOG"
