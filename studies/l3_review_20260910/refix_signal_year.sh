#!/usr/bin/env bash
# Rebuild the signal nominal ROOT files without the `year` branch, then re-hadd and
# re-gate.
#
# Why: `year` is a STRING in the signal parquet. pd.to_numeric(coerce) in the converter
# turned it into NaN/double for 2022-2023 and int64 for 2024, so hadd -- which fixes its
# schema from the FIRST input -- silently dropped the entire 2024 tree from all 14
# signal run3.root files. 23% of the signal statistics, one Warning, exit 0. GATE 4
# caught it; nothing else would have.
#
# Scope: only the 70 signal nominal conversions. bkg/data parquet have no string columns
# and passed GATE 4 untouched. The 1120 signal SYSTEMATIC ROOT files here still carry
# `year`, but they are a dead intermediate -- the MVA-scoring stage reconverts everything
# from parquet into run3_bdt_scored_fsrfix with the fixed converter. Do not hadd them.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
IN=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC
OUT=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix
CONV=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Parque2Root_BDT.py
MASSES="M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M15 M20 M25 M30"
YEARS="2022preEE 2022postEE 2023preBPix 2023postBPix 2024"

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u
export PYTHONPATH="${PYTHONPATH:-}:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA"

echo "=== [1/4] reconverting 70 signal nominal files (5 at a time) ==="
n=0
for m in $MASSES; do
  for y in $YEARS; do
    python $CONV -i "$IN/mA_${m}_${y}/merged_nominal.parquet" \
                 -o "$OUT/mA_${m}/${y}.root" --split \
      > "$D/logs_fsrfix/reconv_mA_${m}_${y}.log" 2>&1 &
    n=$((n+1))
    if [ $((n % 5)) -eq 0 ]; then wait; echo "  ... $n/70 done $(date '+%H:%M:%S')"; fi
  done
done
wait
echo "  ... 70/70 done $(date '+%H:%M:%S')"

echo "=== [2/4] verifying the year branch is gone ==="
python - <<'PY'
import uproot, sys
B="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix"
bad=[]
for m in ["M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M15","M20","M25","M30"]:
    for y in ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]:
        p="%s/mA_%s/%s.root"%(B,m,y)
        try:
            with uproot.open(p) as f:
                if "year" in f["inclusive"]: bad.append(p)
        except Exception as e:
            bad.append("%s (%s)"%(p,str(e)[:50]))
print("files still carrying `year` or unreadable: %d"%len(bad))
for b in bad[:10]: print("   ",b)
sys.exit(1 if bad else 0)
PY
[ $? -ne 0 ] && { echo "STOPPING: year branch survived"; exit 1; }

echo "=== [3/4] re-hadd the signal run3.root files ==="
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile
bash 2_prepare_rootfile.sh --only sig
echo "raw exit: $? (not a verdict)"

echo "=== [4/4] re-running GATE 4 and GATE 4b ==="
python $D/gate4_hadd_closure.py;  g4=$?
python $D/gate4b_branch_closure.py; g4b=$?
echo "GATE 4 exit=$g4   GATE 4b exit=$g4b"
[ $g4 -eq 0 ] && [ $g4b -eq 0 ] && echo "BOTH GATES PASSED" || echo "GATES FAILED"
