#!/usr/bin/env bash
# 等 zee_zmmg control 三個 poller 跑完(pollers gone + queue empty)，再核對每個
# task 的 merged_nominal.parquet 是否齊全，輸出核對報告。
set -uo pipefail
CON='regexp("eos_logs/(mc|data_egamma|data_muon)/", Cmd)'
ROOT=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_zee_zmmg_ctrl
REP=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/merged_verify_zee_zmmg.log
ts(){ date "+%F %T"; }
np(){ ps -u "$USER" -o args 2>/dev/null | grep "[r]un_analysis.py" | grep -c control_zee_zmmg; }
nq(){ condor_q "$USER" -constraint "$CON" -af ClusterId 2>/dev/null | wc -l; }
streak=0; i=0
while [ $i -lt 360 ]; do
  p=$(np); q=$(nq)
  if [ "$p" -eq 0 ] && [ "$q" -eq 0 ]; then streak=$((streak+1)); else streak=0; fi
  [ "$streak" -ge 5 ] && break
  sleep 120; i=$((i+1))
done
{
  echo "===== zee_zmmg merged 核對 $(ts) (pollers=$(np) queue=$(nq)) ====="
  grand_t=0; grand_m=0
  for s in mc data_egamma data_muon; do
    d="$ROOT/$s"
    mapfile -t tasks < <(find "$d" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | grep -E '_(2022|2023|2024)' | sort)
    nt=${#tasks[@]}; nm=0; missing=()
    for t in "${tasks[@]}"; do
      if [ -f "$t/merged_nominal.parquet" ]; then nm=$((nm+1)); else missing+=("$(basename "$t")"); fi
    done
    pj=$(find "$d" -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l)
    echo "--- $s : tasks=$nt merged=$nm perjob_out=$pj"
    [ ${#missing[@]} -gt 0 ] && printf '    MISSING merged: %s\n' "${missing[*]}"
    grand_t=$((grand_t+nt)); grand_m=$((grand_m+nm))
  done
  echo "===== TOTAL tasks=$grand_t merged=$grand_m ====="
  [ "$grand_m" -eq "$grand_t" ] && echo "RESULT: ALL MERGED OK" || echo "RESULT: INCOMPLETE ($((grand_t-grand_m)) task 缺 merged)"
} | tee "$REP"
