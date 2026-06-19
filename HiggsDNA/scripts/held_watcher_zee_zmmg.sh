#!/usr/bin/env bash
# Watch the zee_zmmg control-sample condor jobs; auto-recover OOM (code 34) held
# jobs by escalating RequestMemory and releasing. Logs every action + liveness.
# 記憶體升階：current+4000，上限 24000；超過仍 held → 留著並標記需手動處理。
# 結束條件：poller 都不在了 AND queue 連續 5 次(=10min)看不到 control job 才退出，
#           避免 chunk_size=1 慢速提交時誤判「跑完」。
set -uo pipefail
CON='regexp("eos_logs/(mc|data_egamma|data_muon)/", Cmd)'
ROOT=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_zee_zmmg_ctrl
LOG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/held_watcher_zee_zmmg.log
MAXMEM=24000; STEP=4000
ts(){ date "+%F %T"; }
n_pollers(){ ps -u "$USER" -o args 2>/dev/null | grep "[r]un_analysis.py" | grep -c control_zee_zmmg; }
echo "$(ts) [start] held-watcher v2 for zee_zmmg control" >> "$LOG"
empty_streak=0
i=0
while [ $i -lt 360 ]; do
  nq=$(condor_q "$USER" -constraint "$CON" -af ClusterId 2>/dev/null | wc -l)
  nheld=$(condor_q "$USER" -constraint "$CON && JobStatus==5" -af ClusterId 2>/dev/null | wc -l)
  np=$(n_pollers)
  perjob=$(find "$ROOT" -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l)
  merged=$(find "$ROOT" -name 'merged_nominal.parquet' 2>/dev/null | wc -l)
  echo "$(ts) [tick] pollers=$np inQ=$nq held=$nheld perjob_out=$perjob merged=$merged" >> "$LOG"

  # 處理 OOM held
  for cid in $(condor_q "$USER" -constraint "$CON && JobStatus==5" -af ClusterId 2>/dev/null); do
    code=$(condor_q "$cid" -af HoldReasonCode 2>/dev/null | head -1)
    cur=$(condor_q "$cid" -af RequestMemory 2>/dev/null | head -1)
    task=$(condor_q "$cid" -af Cmd 2>/dev/null | grep -oE "eos_logs/[^/]+/[^/]+" | head -1)
    if [ "$code" = "34" ]; then
      if [ "$cur" -ge "$MAXMEM" ]; then
        echo "$(ts) [MANUAL] $cid ($task) still held at ${cur}MB (>=cap)" >> "$LOG"; continue
      fi
      new=$(( cur + STEP )); [ "$new" -gt "$MAXMEM" ] && new=$MAXMEM
      condor_qedit "$cid" RequestMemory "$new" >/dev/null 2>&1
      condor_release "$cid" >/dev/null 2>&1
      echo "$(ts) [release] $cid ($task) OOM mem ${cur}->${new}MB" >> "$LOG"
    else
      echo "$(ts) [skip] $cid ($task) held code=$code (not OOM) — manual check" >> "$LOG"
    fi
  done

  # 結束判斷：poller 沒了且 queue 空，連續 5 次才退
  if [ "$np" -eq 0 ] && [ "$nq" -eq 0 ]; then
    empty_streak=$((empty_streak+1))
  else
    empty_streak=0
  fi
  if [ "$empty_streak" -ge 5 ]; then
    echo "$(ts) [done] pollers gone & queue empty x5 -> exit (perjob_out=$perjob merged=$merged)" >> "$LOG"; break
  fi
  i=$((i+1)); sleep 120
done
echo "$(ts) [exit] watcher loop ended (iter=$i)" >> "$LOG"
