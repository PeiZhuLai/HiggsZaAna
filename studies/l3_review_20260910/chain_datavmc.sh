#!/usr/bin/env bash
# Stage 6: finish scoring, then run dataVmc (determine MVA cut).
#
# Waits for the two oversized 2024 jobs, gates on GATE 5 over all 1207, hadds the 2024
# DY flavor split, then submits dataVmc on condor.
#
# Why dataVmc must go through condor and not a background shell here: the lxplus login
# node arbiter (CPU quota watchdog) repeatedly kills a full dataVmc reproduction -- the
# first run got lucky and survived 4 h, after which the accumulated badness score meant
# it was killed within 1-5 minutes, and even the sleeping watcher process was killed.
# Heavy work goes to condor; this script only submits and polls.
#
# The env settings in dataVmc.submit were fixed on 2026-09-23 to match what the scoring
# jobs needed (LOCALIZE_CONDA_ENV=1, SETUP_CONDA_ENV=0, CONDA_TARBALL unset) -- the old
# auto/auto plus a nonexistent AFS tarball path is what killed 832 of 1207 scoring jobs.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/Condor

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u

count_bn() { condor_q -constraint "JobBatchName==\"$1\"" -format "%d\n" ClusterId 2>/dev/null | wc -l; }
held_bn()  { condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -format "%d\n" ClusterId 2>/dev/null | wc -l; }

wait_bn() {
  local n held last=-1
  while true; do
    n=$(count_bn "$1"); held=$(held_bn "$1")
    if [ "$n" != "$last" ]; then echo "  [$(date '+%H:%M:%S')] $2: $n in queue ($held held)"; last=$n; fi
    [ "$n" -eq 0 ] && return 0
    if [ "$held" -gt 0 ] && [ "$held" -eq "$n" ]; then
      echo "  [$(date '+%H:%M:%S')] $2: all remaining held"
      condor_q -constraint "JobBatchName==\"$1\" && JobStatus==5" -af HoldReason 2>/dev/null | sort | uniq -c | head -3
      return 1
    fi
    sleep 120
  done
}

echo "=============== [1/5] wait for the two 8 GB 2024 jobs ==============="
wait_bn "p2root_MVAScore_big2" "big2" || { echo "big2 held -- stopping"; exit 1; }

echo
echo "=============== [2/5] GATE 5 over all 1207 scored files ==============="
# since = 2026-09-22 00:00, when run3_bdt_scored_fsrfix was created. Using "today 00:00"
# here previously flagged 352 perfectly good files produced on the 22nd as STALE.
python $D/gate5_scored.py --since $(date -d "2026-09-22 00:00" +%s) \
  || { echo "GATE 5 FAILED -- not starting dataVmc"; exit 1; }

echo
echo "=============== [3/5] hadd the 2024 DY flavor split ==============="
cd "$C"
bash 4_prepaare_2024DYJetsToLL.sh
echo "raw exit: $? (not a verdict)"
python - <<'PY'
import uproot, os, sys
B="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
out="%s/DYJetsToLL/2024.root" % B
ins=["%s/DYJetsTo2%s/2024.root" % (B,c) for c in ("E","Mu","Tau")]
if not os.path.exists(out):
    print("MISSING", out); sys.exit(1)
o=uproot.open(out)["inclusive"].num_entries
t=0
for p in ins:
    if not os.path.exists(p): print("MISSING INPUT", p); sys.exit(1)
    t+=uproot.open(p)["inclusive"].num_entries
print("hadd closure: output=%d  sum(inputs)=%d  %s" % (o,t,"OK" if o==t else "LOST %+d"%(o-t)))
sys.exit(0 if o==t else 1)
PY
[ $? -ne 0 ] && { echo "2024 DY hadd closure FAILED -- stopping"; exit 1; }

echo
echo "=============== [4/5] submit dataVmc ==============="
cd "$P"
SKIP_ENV_PACK=0 NO_SUBMIT=1 bash 1_submit_dataVmc_condor.sh
rows=$(grep -cve '^\s*$' dataVmc_jobs.txt 2>/dev/null || echo 0)
echo "  dataVmc_jobs.txt rows: $rows"
[ "$rows" -eq 0 ] && { echo "no dataVmc jobs generated -- stopping"; exit 1; }
condor_submit dataVmc.submit
wait_bn "dataVmc_fsrfix" "dataVmc" || echo "  (some held -- merge step will show what landed)"

echo
echo "=============== [5/5] merge + draw ==============="
bash 3_merge_dataVmc_condor.sh
echo "raw exit: $? (not a verdict)"
echo "=============== DATAVMC CHAIN DONE $(date '+%F %T') ==============="
