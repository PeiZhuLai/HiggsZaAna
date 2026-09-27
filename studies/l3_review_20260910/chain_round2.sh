#!/usr/bin/env bash
# Round 2 of the FSR-fix chain: rerun p2root with (a) the fixed converter and (b) the
# regenerated sideband reweight JSON, then re-hadd and re-gate, then check that the
# sideband reweight has converged.
#
# Why a second round:
#  1. `year` -- a string column in the signal parquet -- became NaN/double for 2022-2023
#     and int64 for 2024, so hadd dropped the whole 2024 tree from every signal
#     run3.root (-23%). The converter now drops it, but the 1120 signal SYSTEMATIC ROOT
#     files from round 1 still carry it. Rerunning cleans all 1207 in one pass.
#  2. The sideband JSON is derived from p2root output, yet p2root stamps
#     weight_sideband_rwgt into every ROOT file from that same JSON. Round 1 used the
#     2026-06-03 JSON (pre-FSR-fix). The trainer reads the JSON while dataVmc reads the
#     ROOT branch, so the two disagreed. Round 2 makes them the same object.
#     Shape factors moved by rms 0.26-1.22% between the two JSONs, so round 2 is
#     expected to CONVERGE, not to move again -- that is what step 6 verifies.
#
# Fail-closed: every gate must pass before the next step runs.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
M=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA
START=$(date +%s)

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"
echo "[env] python=$(which python)"
echo "[env] production sideband JSON: $(python -c "import json;print(json.load(open('$M/reweights/sideband_run3_iterative.json'))['created'])")"

echo
echo "=============== [1/6] p2root, all 1207 conversions ==============="
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
bash 1_run_P2Root.sh
echo "raw exit: $? (not a verdict)"

echo
echo "=============== [2/6] GATE 3: merged rows == inclusive entries ==============="
python $D/gate3_reconcile.py --since $START || { echo "GATE 3 FAILED -- stopping"; exit 1; }

echo
echo "=============== [3/6] verify the year branch is gone everywhere ==============="
python - <<'PY'
import uproot, glob, sys, os
B="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
bad=[]
files=[p for p in glob.glob(B+"/*/*.root") if os.path.basename(p)!="run3.root"]
for p in files:
    try:
        with uproot.open(p) as f:
            if "year" in f["inclusive"]: bad.append(p)
    except Exception as e:
        bad.append("%s (%s)"%(p,str(e)[:50]))
print("checked %d files; still carrying `year` or unreadable: %d" % (len(files), len(bad)))
for b in bad[:10]: print("   ", b)
sys.exit(1 if bad else 0)
PY
[ $? -ne 0 ] && { echo "STOPPING: year branch survived"; exit 1; }

echo
echo "=============== [4/6] hadd ==============="
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
bash 2_prepare_rootfile.sh
echo "raw exit: $? (not a verdict)"

echo
echo "=============== [5/6] GATE 4 + GATE 4b ==============="
python $D/gate4_hadd_closure.py   || { echo "GATE 4 FAILED -- stopping";  exit 1; }
python $D/gate4b_branch_closure.py || { echo "GATE 4b FAILED -- stopping"; exit 1; }

echo
echo "=============== [6/6] sideband convergence check ==============="
cd $M/scripts
python3 1_make_sideband_reweight.py --output-json $M/reweights/sideband_run3_iterative_round2.json || exit 1
echo
echo "--- round 1 JSON (now in production) vs round 2 ---"
echo "--- these should agree to well under the round-0 -> round-1 move ---"
python3 $D/compare_sideband_json.py \
  $M/reweights/sideband_run3_iterative.json \
  $M/reweights/sideband_run3_iterative_round2.json

echo
echo "=============== ROUND 2 COMPLETE $(date '+%F %T') ==============="
