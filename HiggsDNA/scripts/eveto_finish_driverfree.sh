#!/usr/bin/env bash
# HZa photon e-veto 研究 收尾（不靠 HiggsDNA driver，避開 driver hang）：
#   等我的 eveto condor jobs drain → 掃 NPS 缺 nominal parquet 的 job → 直接 condor_submit
#   其既有 batch_submit 檔補投 → 迴圈到全齊 → collect cutflow → plot。
# 只碰 eos_logs/{Sig_MC,Bkg_MC} 的 pelai jobs（不動 user 其他 production）。
set -uo pipefail

REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
PLOTREPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
EOSLOG="$REPO/eos_logs"
NPS=/eos/cms/store/group/phys_susy/NPS-25-014/HZa/parquet_eveto
LOGDIR="$REPO/eveto_submit_logs"
STAMP=$(date +%Y%m%d_%H%M%S)
FLOG="$LOGDIR/finish_${STAMP}.log"
CONDA_SH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
MAXPASS=12
mkdir -p "$LOGDIR"
log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$FLOG"; }

active_cnt(){ condor_q pelai -af JobStatus Cmd 2>/dev/null \
  | grep -E "eos_logs/(Sig_MC|Bkg_MC)/" | awk '$1==1||$1==2' | wc -l; }

wait_drain(){
  local st=0 a
  while :; do
    a=$(active_cnt)
    log "  active(idle+run)=$a"
    if [ "$a" -eq 0 ]; then st=$((st+1)); [ "$st" -ge 3 ] && break; else st=0; fi
    sleep 300
  done
}

# 列出缺 nominal parquet 的 job 的 batch_submit 檔（missing = 失敗/未跑）
list_missing_submits(){  # -> stdout: 每行一個 batch_submit.txt 路徑
  local leaf sampdir jobdir jname sub npsjob
  for leaf in Sig_MC Bkg_MC; do
    for jobdir in "$EOSLOG/$leaf"/*/job_*; do
      [ -d "$jobdir" ] || continue
      jname=$(basename "$jobdir")                 # job_N
      sampdir=$(basename "$(dirname "$jobdir")")  # <sample>_<year>
      npsjob="$NPS/$leaf/$sampdir/$jname"
      # done 判準：NPS job dir 有 *nominal.parquet
      if ! ls "$npsjob"/*nominal.parquet >/dev/null 2>&1; then
        sub=$(ls "$jobdir"/*batch_submit*.txt 2>/dev/null | head -1)
        [ -n "$sub" ] && echo "$sub"
      fi
    done
  done
}

log "driver-free finisher start (stamp=$STAMP)."
log "waiting for current eveto batch to drain..."
sleep 30
wait_drain
log "current batch drained."

for PASS in $(seq 1 $MAXPASS); do
  log "=== pass $PASS: scanning NPS for missing nominal parquet ==="
  mapfile -t MISS < <(list_missing_submits)
  n=${#MISS[@]}
  log "pass $PASS: missing jobs = $n"
  if [ "$n" -eq 0 ]; then log "all jobs complete."; break; fi
  # 直接 condor_submit 每個缺件 job 的既有 submit 檔
  ok=0
  for sub in "${MISS[@]}"; do
    if condor_submit "$sub" >>"$FLOG" 2>&1; then ok=$((ok+1)); fi
  done
  log "pass $PASS: resubmitted $ok / $n jobs. waiting to drain..."
  sleep 60
  wait_drain
done

# ---- collect ----
log "running collector (hza_ana) ..."
set +u; source "$CONDA_SH"; conda activate hza_ana; set -u
export PYTHONPATH="$REPO"; cd "$REPO"
python scripts/5_collect_cutflow.py >>"$FLOG" 2>&1
python - >>"$FLOG" 2>&1 <<'PY'
import json,glob
fs=sorted(glob.glob("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list/cutflow_Sig_MC_*_2022preEE.json"))
if fs:
    d=json.load(open(fs[0])); cf=d.get("cutflows",{})
    for k in ("zgammas_ph_eveto_w","zgammas_ph_noeveto_w"):
        c=cf.get(k,{}); print("  ",fs[0].split('/')[-1],k,"allcuts=",c.get("all cuts") or c.get("event"))
PY
set +u; conda deactivate; set -u

# ---- plot ----
log "running plotter (higgs-zg-plot_python) ..."
set +u; source "$CONDA_SH"; conda activate higgs-zg-plot_python; set -u
cd "$PLOTREPO"
python scripts/plot_evetosigniVmA.py >>"$FLOG" 2>&1
log "plotter exit=$?"
ls -la "$PLOTREPO/plots/evetosigniVmA/" >>"$FLOG" 2>&1
log "finisher DONE. log=$FLOG"
