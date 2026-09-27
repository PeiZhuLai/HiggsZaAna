#!/usr/bin/env bash
# Stage 5: BDT scoring on condor. Smoketest 4 jobs first, gate on the PRODUCT, then
# submit the full 1207.
#
# Two failure modes this guards against, neither visible in an exit code:
#  * add_mva_scores() catches a model-load failure, fills every score with NaN, prints a
#    WARNING and returns. The job exits 0 and writes a valid ROOT file with no usable
#    scores. GATE 5 reads the score branches and rejects all-NaN.
#  * The conda env is fetched per job by xrdcp from EOS (the packer deletes the local
#    copy after upload). If that path is wrong every job dies before running anything.
#    The smoketest is what proves the fetch works on a worker, not on lxplus.
#
# Batch names are exact -- the schedd is shared, and a loose constraint has previously
# matched another session's jobs and made a wait loop never finish.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor
TEST_BN="p2root_MVAScore_TEST"
FULL_BN="p2root_MVAScore"

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u

count_jobs() { condor_q -constraint "JobBatchName==\"$1\"" -format "%d\n" ClusterId 2>/dev/null | wc -l; }
count_held() { condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -format "%d\n" ClusterId 2>/dev/null | wc -l; }

wait_for() {   # $1 = batch name, $2 = label
  local n held last=-1
  while true; do
    n=$(count_jobs "$1"); held=$(count_held "$1")
    if [ "$n" != "$last" ]; then
      echo "  [$(date '+%H:%M:%S')] $2: $n in queue ($held held)"; last=$n
    fi
    [ "$n" -eq 0 ] && break
    if [ "$held" -gt 0 ] && [ "$held" -eq "$n" ]; then
      echo "  [$(date '+%H:%M:%S')] $2: ALL remaining jobs held -- stopping the wait"
      condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -af HoldReason 2>/dev/null | sort | uniq -c | head
      return 1
    fi
    sleep 60
  done
  return 0
}

echo "=============== [1/5] preflight: three models must be fresh ==============="
python - <<'PY'
import os, pickle, sys, time
U="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using"
cut = time.time() - 24*3600
exp = {"model_Za_BDT_run3.pkl":16, "model_Za_BDT_lowmass_run3.pkl":5, "model_Za_BDT_highmass_run3.pkl":16}
ok=True
for n,nf_exp in exp.items():
    p=os.path.join(U,n)
    if not os.path.exists(p): print("MISSING",p); ok=False; continue
    mt=os.path.getmtime(p)
    m=pickle.load(open(p,"rb")); nf=getattr(m,"n_features_in_",None)
    good = mt>=cut and nf==nf_exp
    print("%-34s %-4s mtime=%s n_features=%s" % (n,"OK" if good else "BAD",
          time.strftime('%F %T',time.localtime(mt)), nf))
    ok &= good
sys.exit(0 if ok else 1)
PY
[ $? -ne 0 ] && { echo "STOPPING: models are not all freshly trained"; exit 1; }

cd "$C"
echo
echo "=============== [2/5] smoketest: 4 representative jobs ==============="
echo "  $(wc -l < joblist_test.tsv) jobs: signal nominal, signal syst, bkg, data"
SMOKE_START=$(date +%s)
condor_submit 2_submit_test.sub
wait_for "$TEST_BN" "smoketest" || { echo "smoketest jobs held -- stopping"; exit 1; }

echo
echo "=============== [3/5] GATE 5 + mechanism check on the smoketest ==============="
python $D/gate5_scored.py --joblist $C/joblist_test.tsv --since $SMOKE_START \
  || { echo "SMOKETEST FAILED GATE 5 -- not submitting the full set"; exit 1; }

# A 4-job smoketest passing is NOT evidence that the full set will pass. On 2026-09-22
# the smoketest passed 4/4 and then 832 of 1207 jobs died with EOS "Input/output error",
# because with SETUP_CONDA_ENV=auto each job activated the conda env directly off EOS
# before localizing. Four jobs do not overload EOS; twelve hundred do. So verify the
# MECHANISM in the job logs, not just the products: the env must come from the xrdcp'd
# tarball and the job must never touch the EOS env.
echo "--- mechanism: conda env must be localized, EOS env must not be activated ---"
smoke_cluster=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$(tail -20 $D/logs_fsrfix/scoring.log)" | tail -1)
nfetch=$(grep -l "Fetch conda tarball" $C/logs/${smoke_cluster}.*.out 2>/dev/null | wc -l)
nact=$(grep -l "Activate conda env"  $C/logs/${smoke_cluster}.*.out 2>/dev/null | wc -l)
njobs=$(wc -l < $C/joblist_test.tsv)
echo "  cluster=$smoke_cluster  jobs=$njobs  localized=$nfetch  activated-from-EOS=$nact"
if [ "$nfetch" -ne "$njobs" ] || [ "$nact" -ne 0 ]; then
  echo "  MECHANISM CHECK FAILED: every job must localize ($nfetch/$njobs) and none may"
  echo "  activate from EOS ($nact must be 0). Not submitting the full set."
  exit 1
fi
echo "  mechanism OK"

echo
echo "=============== [4/5] full submission: 1207 jobs ==============="
FULL_START=$(date +%s)
echo "  joblist lines: $(wc -l < joblist.tsv)"
condor_submit 2_submit.sub
wait_for "$FULL_BN" "full" || echo "  (wait ended early -- GATE 5 below will show what landed)"

echo
echo "=============== [5/5] GATE 5 on the full set ==============="
# since = today 00:00, not FULL_START: run3_bdt_scored_fsrfix is a directory created
# today, so there are no stale products to guard against, and files written by an earlier
# batch this same day are legitimately good.
python $D/gate5_scored.py --joblist $C/joblist.tsv --since $(date -d "today 00:00" +%s)
rc=$?
echo "GATE 5 exit=$rc"
echo "=============== SCORING CHAIN DONE $(date '+%F %T') ==============="
exit $rc
