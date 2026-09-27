#!/usr/bin/env bash
# Stage 5b: resubmit the 840 scoring jobs that died on EOS I/O, with the env fix.
#
# What went wrong the first time: SETUP_CONDA_ENV=auto made every job activate the conda
# env directly off EOS, and only THEN decide to localize -- the "auto" localize test is
# whether CONDA_PREFIX starts with /eos, which it cannot know before activating from
# /eos. 1207 jobs doing that at once gave 832 x "Input/output error" (Errno 5),
# return value 1, ZERO held, and a queue that drained normally. Nothing looked wrong
# except that 840 ROOT files were never written.
#
# The 4-job smoketest passed 4/4 immediately before that. Four jobs do not overload EOS.
# That is why step [3] now checks the MECHANISM in the job logs -- every job must fetch
# the tarball, and none may activate from EOS -- instead of only checking products.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor
TEST_BN="p2root_MVAScore_TEST"
RESUB_BN="p2root_MVAScore_resub"

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u

count_jobs() { condor_q -constraint "JobBatchName==\"$1\"" -format "%d\n" ClusterId 2>/dev/null | wc -l; }
count_held() { condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -format "%d\n" ClusterId 2>/dev/null | wc -l; }

wait_for() {
  local n held last=-1
  while true; do
    n=$(count_jobs "$1"); held=$(count_held "$1")
    if [ "$n" != "$last" ]; then echo "  [$(date '+%H:%M:%S')] $2: $n in queue ($held held)"; last=$n; fi
    [ "$n" -eq 0 ] && break
    if [ "$held" -gt 0 ] && [ "$held" -eq "$n" ]; then
      echo "  [$(date '+%H:%M:%S')] $2: ALL remaining held -- stopping"
      condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -af HoldReason 2>/dev/null | sort | uniq -c | head
      return 1
    fi
    sleep 60
  done
  return 0
}

check_mechanism() {   # $1 = cluster, $2 = expected job count, $3 = label
  local nfetch nact
  nfetch=$(grep -l "Fetch conda tarball" $C/logs/$1.*.out 2>/dev/null | wc -l)
  nact=$(grep -l "Activate conda env"  $C/logs/$1.*.out 2>/dev/null | wc -l)
  echo "  $3: cluster=$1 jobs=$2 localized=$nfetch activated-from-EOS=$nact"
  if [ "$nfetch" -ne "$2" ] || [ "$nact" -ne 0 ]; then
    echo "  MECHANISM CHECK FAILED (need localized=$2, activated=0)"
    return 1
  fi
  echo "  mechanism OK"
  return 0
}

cd "$C"
echo "=============== [1/4] smoketest with the env fix ==============="
NTEST=$(wc -l < joblist_test.tsv)
SMOKE_START=$(date +%s)
out=$(condor_submit 2_submit_test.sub 2>&1); echo "$out"
SMOKE_CLUSTER=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out" | tail -1)
wait_for "$TEST_BN" "smoketest" || { echo "smoketest held -- stopping"; exit 1; }

echo
echo "=============== [2/4] smoketest gates: product AND mechanism ==============="
python $D/gate5_scored.py --joblist $C/joblist_test.tsv --since $SMOKE_START \
  || { echo "SMOKETEST FAILED GATE 5 -- not resubmitting"; exit 1; }
check_mechanism "$SMOKE_CLUSTER" "$NTEST" "smoketest" || exit 1

echo
echo "=============== [3/4] resubmit the 840 failed jobs ==============="
NRESUB=$(wc -l < joblist_resub.tsv)
echo "  joblist_resub.tsv: $NRESUB jobs, max_materialize=150"
out=$(condor_submit 2_submit_resub.sub 2>&1); echo "$out"
RESUB_CLUSTER=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out" | tail -1)
wait_for "$RESUB_BN" "resub" || echo "  (wait ended early -- GATE 5 below shows what landed)"
check_mechanism "$RESUB_CLUSTER" "$NRESUB" "resub" || echo "  (mechanism imperfect -- see GATE 5)"

echo
echo "=============== [4/4] GATE 5 on all 1207 ==============="
python $D/gate5_scored.py --joblist $C/joblist.tsv --since $(date -d "today 00:00" +%s)
rc=$?
echo "GATE 5 exit=$rc"
echo "=============== RESUB CHAIN DONE $(date '+%F %T') ==============="
exit $rc
