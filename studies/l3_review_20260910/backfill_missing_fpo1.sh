#!/usr/bin/env bash
#
# Resubmit the fpo=1 jobs that produced no chunk.
#
# WHY NOT THE DRIVER
#   run_analysis.py retires a job after 5 worker failures and then counts it as
#   complete -- that is how the stage reported "100 percent" while 2400 chunks
#   were missing. UNRETIRE_JOBS=1 would revive them, but the driver also
#   re-walks 28k job configs on EOS (that took ~45 min last time) before it
#   submits anything. The per-job .sub files already exist, already carry
#   fpo=1 / 3000 MB / 1 cpu / workday, and PYTHONPATH picks up the current code,
#   so submitting them directly is both faster and narrower.
#
# WHY THE EXPORTS
#   submit_template.txt is getenv = True, so a job inherits the SUBMITTING
#   shell's environment. Submitting these files from a bare shell gave 20/20
#   jobs exit 127 in 3 s on 2026-09-15: the wrapper leaves "python" literal and
#   without HZA_BMODE_PACK nothing ever activates an environment.
#
# INPUT   logs_fsrfix/missing_jobs.txt   ("<stage> <sample> <job_N>" per line,
#                                         written by gap_report_fpo1.py)
#
# USAGE
#   bash backfill_missing_fpo1.sh             # dry run: counts and one example
#   DRY_RUN=0 bash backfill_missing_fpo1.sh   # submit
#   DRY_RUN=0 MAX=200 bash ...                # submit only the first 200
#
set -uo pipefail
DRY_RUN="${DRY_RUN:-1}"
MAX="${MAX:-0}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
R=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
# Default to the v2 list: the v1 list counted parquets, and at fpo=1 "no parquet"
# usually means zero events passed selection, not a failed job. v1 called 2286
# healthy jobs missing and caused 200 needless resubmissions. The completion
# marker is the summary json.
LIST="${LIST:-$D/logs_fsrfix/missing_jobs_v2.txt}"
LOG=$D/logs_fsrfix/backfill_missing.log
ENVDIR=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana

[ -s "$LIST" ] || { echo "[ERROR] $LIST is empty or absent -- run gap_report_fpo1.py first"; exit 1; }

export CONDA_PREFIX="$ENVDIR"
export PATH="$ENVDIR/bin:$PATH"
export PYTHONPATH="$R"
export X509_USER_PROXY="${X509_USER_PROXY:-/tmp/x509up_u$(id -u)}"
export HZA_BMODE_PACK="${HZA_BMODE_PACK:-root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz}"
export HZA_STAGE_INPUTS=1

left=$(voms-proxy-info -timeleft 2>/dev/null || echo 0)
echo "proxy left: ${left}s"
[ "${left:-0}" -lt 86400 ] && { echo "[ERROR] proxy under 24 h -- refresh it AND copy it to the repo"; exit 1; }
cp -p "$X509_USER_PROXY" "$R/x509up_u$(id -u)" 2>/dev/null
"$ENVDIR/bin/openssl" x509 -noout -enddate -in "$R/x509up_u$(id -u)"

n=$(wc -l < "$LIST"); echo "missing jobs listed: $n"
subs=(); nosub=0
while read -r stage smp job; do
    [ -z "${job:-}" ] && continue
    f="$R/eos_logs/$stage/$smp/$job/${smp}_batch_submit_${job/job_/job}.txt"
    if [ -f "$f" ]; then subs+=("$f"); else nosub=$((nosub+1)); fi
done < "$LIST"
echo "submit files found: ${#subs[@]}   missing a .sub: ${nosub}"
[ "${#subs[@]}" -eq 0 ] && { echo "[ERROR] no submit files resolved -- check the path pattern"; exit 1; }
echo "example: ${subs[0]}"
grep -hE 'RequestMemory|RequestCpus|JobFlavour' "${subs[0]}" | sed 's/^/  /'

[ "$MAX" -gt 0 ] && subs=("${subs[@]:0:$MAX}")
if [ "$DRY_RUN" = "1" ]; then
    echo "would condor_submit ${#subs[@]} jobs"; exit 0
fi
ok=0; bad=0
for f in "${subs[@]}"; do
    if condor_submit "$f" >/dev/null 2>&1; then ok=$((ok+1)); else bad=$((bad+1)); fi
done
echo "$(date '+%F %T') backfill submitted ok=${ok} failed=${bad} of ${#subs[@]}" | tee -a "$LOG"
