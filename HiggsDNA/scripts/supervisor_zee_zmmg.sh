#!/usr/bin/env bash
# Auto-supervisor for the zee_zmmg control-sample production.
# 反覆卡 framework monitoring loop -> 自動 kill+resume 重投未完成 job、處理 OOM hold、
# 把「該 task 全部 job 都有 output」的 task 用 pyarrow 自行 merge（繞過會 hang 的框架 merge）。
# 全部 task 都 merged 或達 cycle 上限才結束。
set -uo pipefail
export CONDA_PREFIX=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
export PATH=$CONDA_PREFIX/bin:$PATH
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
ROOT=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_zee_zmmg_ctrl
LOG=$REPO/eos_logs/supervisor_zee_zmmg.log
cd "$REPO"

CYCLE=300; MAXCYCLE=320; STALE_MIN=20; MEMCAP=24000; MEMSTEP=4000; RESTCAP=12
YEARS=2022preEE,2022postEE,2023preBPix,2023postBPix,2024
declare -A SCRIPT=( [mc]=run_ana_control_zee_zmmg.sh [data_egamma]=run_ana_data_control_zee_zmmg.sh [data_muon]=run_ana_data_control_zee_zmmg.sh )
declare -A SLIST=( [mc]="" [data_egamma]=Data [data_muon]=Data_tnp_zmmg )
declare -A REST=( [mc]=0 [data_egamma]=0 [data_muon]=0 )
declare -A LASTDONE=( [mc]=-1 [data_egamma]=-1 [data_muon]=-1 )
declare -A GIVEUP=( [mc]=0 [data_egamma]=0 [data_muon]=0 )
declare -A TLASTGOT=(); declare -A TSTALL=()   # per-task plateau tracking
PLATEAU=4        # cycles of no per-task progress before allowing partial merge
PARTIAL_PCT=97   # merge partial only if got/exp >= this %

ts(){ date "+%F %T"; }
say(){ echo "$(ts) $*" >> "$LOG"; }
conq(){ condor_q "$USER" -constraint "regexp(\"eos_logs/$1/\", Cmd)" -af ClusterId 2>/dev/null | wc -l; }
poller_pids(){ ps -u "$USER" -o pid,args 2>/dev/null | grep "[r]un_analysis.py" | grep "$ROOT/$1 " | awk '{print $1}'; }
total_jobs(){ grep -hoE "split over [0-9]+ total jobs" "$REPO/eos_logs/resub3_$1.log" 2>/dev/null | grep -oE "[0-9]+" | tail -1; }
done_jobs(){ find "$ROOT/$1" -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l; }
log_age_min(){ local f="$REPO/eos_logs/resub3_$1.log"; [ -f "$f" ] || { echo 99999; return; }; echo $(( ( $(date +%s) - $(stat -c %Y "$f") ) / 60 )); }

handle_oom(){
  local reg=$1
  for cid in $(condor_q "$USER" -constraint "regexp(\"eos_logs/$reg/\", Cmd) && JobStatus==5" -af ClusterId 2>/dev/null); do
    local code cur new task
    code=$(condor_q "$cid" -af HoldReasonCode 2>/dev/null | head -1)
    cur=$(condor_q "$cid" -af RequestMemory 2>/dev/null | head -1)
    task=$(condor_q "$cid" -af Cmd 2>/dev/null | grep -oE "eos_logs/[^/]+/[^/]+" | head -1)
    if [ "$code" = "34" ]; then
      if [ "${cur:-0}" -ge "$MEMCAP" ]; then say "[OOM-MANUAL] $cid ($task) stuck at ${cur}MB"; continue; fi
      new=$(( cur + MEMSTEP )); [ "$new" -gt "$MEMCAP" ] && new=$MEMCAP
      condor_qedit "$cid" RequestMemory "$new" >/dev/null 2>&1; condor_release "$cid" >/dev/null 2>&1
      say "[OOM-release] $cid ($task) ${cur}->${new}MB"
    else
      say "[held-other] $cid ($task) code=$code (manual)"
    fi
  done
}

merge_task(){ # $1=region $2=taskdir
  local reg=$1 td=$2
  "$CONDA_PREFIX/bin/python" - "$td" >>"$LOG" 2>&1 <<'PY'
import sys, glob, os, pyarrow.parquet as pq
td=sys.argv[1]
files=sorted(glob.glob(td+'/**/output_job_*_nominal.parquet', recursive=True))
if not files:
    print("  [merge] no per-job files in", td); sys.exit(2)
out=td+'/merged_nominal.parquet'; tmp=out+'.tmp'
w=None; n=0
for f in files:
    try:
        t=pq.read_table(f)
    except Exception as e:
        print("  [merge] skip bad", f, e); continue
    if w is None: w=pq.ParquetWriter(tmp, t.schema)
    try: w.write_table(t); n+=1
    except Exception as e: print("  [merge] schema-skip", f, e)
if w: w.close(); os.replace(tmp, out); print(f"  [merge] {os.path.basename(td)}: {n}/{len(files)} -> merged_nominal.parquet")
else: print("  [merge] nothing merged for", td)
PY
}

restart_region(){
  local reg=$1
  for p in $(poller_pids "$reg"); do kill "$p" 2>/dev/null; done
  sleep 3
  local extra=""; [ -n "${SLIST[$reg]}" ] && extra="SAMPLE_LIST=${SLIST[$reg]}"
  nohup env CLEAN_ANALYSIS_STATE=0 UNRETIRE_JOBS=1 RECONFIGURE_JOBS=0 MERGE_OUTPUTS=0 \
     CONDOR_SUBMIT_CHUNK_SIZE=500 BATCH_SYSTEM=condor YEARS="$YEARS" \
     OUTDIR="$ROOT/$reg" $extra \
     bash "scripts/${SCRIPT[$reg]}" > "$REPO/eos_logs/resub3_$reg.log" 2>&1 &
  REST[$reg]=$(( ${REST[$reg]} + 1 ))
  say "[restart] $reg (attempt ${REST[$reg]})"
}

say "===== supervisor start ====="
cyc=0
while [ $cyc -lt $MAXCYCLE ]; do
  alldone=1
  for reg in mc data_egamma data_muon; do
    [ "${GIVEUP[$reg]}" = "1" ] && continue
    tot=$(total_jobs "$reg"); tot=${tot:-0}
    dn=$(done_jobs "$reg")
    q=$(conq "$reg"); np=$(poller_pids "$reg" | wc -l); age=$(log_age_min "$reg")
    say "[tick] $reg done=$dn/$tot inQ=$q pollers=$np logAge=${age}m restarts=${REST[$reg]}"

    handle_oom "$reg"

    # expected tasks come from submission log dirs (all tasks, even未產出者)
    mapfile -t TASKS < <(find "$REPO/eos_logs/$reg" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | grep -E '_(2022|2023|2024)' | sort)
    ntask=${#TASKS[@]}; nmerg=0
    for ltd in "${TASKS[@]}"; do
      base=$(basename "$ltd"); td="$ROOT/$reg/$base"; key="$reg/$base"
      if [ -f "$td/merged_nominal.parquet" ]; then nmerg=$((nmerg+1)); continue; fi
      exp=$(find "$ltd" -mindepth 1 -maxdepth 1 -type d -name 'job_*' 2>/dev/null | wc -l)
      got=$(find "$td" -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l)
      [ "$exp" -le 0 ] && continue
      # plateau tracking
      if [ "${TLASTGOT[$key]:--1}" -eq "$got" ]; then TSTALL[$key]=$(( ${TSTALL[$key]:-0} + 1 )); else TSTALL[$key]=0; fi
      TLASTGOT[$key]=$got
      if [ "$got" -ge "$exp" ]; then
        merge_task "$reg" "$td"; [ -f "$td/merged_nominal.parquet" ] && nmerg=$((nmerg+1))
      elif [ "${TSTALL[$key]:-0}" -ge "$PLATEAU" ] && [ $(( got*100 )) -ge $(( exp*PARTIAL_PCT )) ]; then
        say "[partial-merge] $key got=$got/$exp (stalled ${TSTALL[$key]} cyc, >=${PARTIAL_PCT}%)"
        merge_task "$reg" "$td"; [ -f "$td/merged_nominal.parquet" ] && nmerg=$((nmerg+1))
      fi
    done
    if [ "$ntask" -gt 0 ] && [ "$nmerg" -ge "$ntask" ]; then
      say "[region-done] $reg merged=$nmerg/$ntask"
      for p in $(poller_pids "$reg"); do kill "$p" 2>/dev/null; done
      continue
    fi
    alldone=0

    # stuck? restart
    if [ "$np" -eq 0 ] || { [ "$age" -ge "$STALE_MIN" ] && [ "$q" -eq 0 ]; }; then
      if [ "$dn" -le "${LASTDONE[$reg]}" ] && [ "${REST[$reg]}" -ge "$RESTCAP" ]; then
        say "[giveup] $reg no progress after ${REST[$reg]} restarts (done=$dn)"; GIVEUP[$reg]=1
      else
        restart_region "$reg"
      fi
    fi
    LASTDONE[$reg]=$dn
  done

  [ "$alldone" = "1" ] && { say "[ALL-DONE] every region merged"; break; }
  cyc=$((cyc+1)); sleep $CYCLE
done

# final report
{
  echo "===== FINAL $(ts) ====="
  gt=0; gm=0
  for reg in mc data_egamma data_muon; do
    mapfile -t FT < <(find "$REPO/eos_logs/$reg" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | grep -E '_(2022|2023|2024)' | sort)
    nt=${#FT[@]}; nm=0; miss=""
    for ltd in "${FT[@]}"; do base=$(basename "$ltd"); if [ -f "$ROOT/$reg/$base/merged_nominal.parquet" ]; then nm=$((nm+1)); else miss="$miss $base"; fi; done
    pj=$(find "$ROOT/$reg" -name 'output_job_*_nominal.parquet' 2>/dev/null | wc -l)
    echo "--- $reg: merged=$nm/$nt perjob=$pj giveup=${GIVEUP[$reg]} restarts=${REST[$reg]}"
    [ -n "$miss" ] && echo "    MISSING merged:$miss"
    gt=$((gt+nt)); gm=$((gm+nm))
  done
  echo "===== TOTAL merged=$gm/$gt ====="
  [ "$gm" -ge "$gt" ] && [ "$gt" -gt 0 ] && echo "RESULT: ALL MERGED OK" || echo "RESULT: INCOMPLETE"
} | tee -a "$LOG"
say "===== supervisor exit (cyc=$cyc) ====="
