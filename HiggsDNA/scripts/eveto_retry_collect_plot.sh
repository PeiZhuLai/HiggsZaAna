#!/usr/bin/env bash
# HZa photon e-veto 研究 orchestrator：
#   等現有批次 drain → 用 --unretire 重投缺 output 的失敗 job(EOS conda transient)
#   → 迴圈直到全齊(或 max passes) → collect cutflow → plot。
# 只碰 eos_logs/Sig_MC 與 eos_logs/Bkg_MC 的 pelai jobs（不動 user 其他 job）。
set -uo pipefail

REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
PLOTREPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
LOGDIR="$REPO/eveto_submit_logs"
STAMP=$(date +%Y%m%d_%H%M%S)
OLOG="$LOGDIR/orch_${STAMP}.log"
mkdir -p "$LOGDIR"
log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$OLOG"; }

CONDA_SH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
HZAENV=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
NPS=/eos/cms/store/group/phys_susy/NPS-25-014/HZa/parquet_eveto
MAXPASS=10

active_cnt(){ condor_q pelai -af JobStatus Cmd 2>/dev/null \
  | grep -E "eos_logs/(Sig_MC|Bkg_MC)/" | awk '$1==1||$1==2' | wc -l; }
inqueue_cnt(){ condor_q pelai -af Cmd 2>/dev/null | grep -cE "eos_logs/(Sig_MC|Bkg_MC)/"; }

wait_drain(){
  local st=0
  while :; do
    local a; a=$(active_cnt)
    log "  active(idle+run)=$a held=$(condor_q pelai -af JobStatus Cmd 2>/dev/null | grep -E 'eos_logs/(Sig_MC|Bkg_MC)/' | awk '$1==5' | wc -l)"
    if [ "$a" -eq 0 ]; then st=$((st+1)); [ "$st" -ge 3 ] && break; else st=0; fi
    sleep 300
  done
}

resubmit(){  # $1=signal|bkgmc  returns 0 if it submitted >0 jobs
  local which="$1" leaf log0
  [ "$which" = signal ] && leaf=Sig_MC || leaf=Bkg_MC
  log0="$LOGDIR/resubmit_${which}_${STAMP}_p${PASS}.log"
  ( cd "$REPO"
    export CONDA_PREFIX="$HZAENV"; export PATH="$HZAENV/bin:$PATH"
    export PYTHONPATH="$REPO"
    CLEAN_ANALYSIS_STATE=0 UNRETIRE_JOBS=1 RETIRE_JOBS=0 MERGE_OUTPUTS=0 \
      BATCH_SYSTEM=condor OUTDIR="$NPS/$leaf" \
      bash scripts/run_ana_${which}.sh --no_systematics >"$log0" 2>&1 )
  local n; n=$(grep -c "Submitted .* jobs in" "$log0" 2>/dev/null)
  log "  resubmit $which: submitted-chunks=$n (log=$log0)"
  [ "${n:-0}" -gt 0 ]
}

log "orchestrator start (stamp=$STAMP). waiting for current batch to drain..."
sleep 60
wait_drain
log "current batch drained."

for PASS in $(seq 1 $MAXPASS); do
  log "=== retry pass $PASS ==="
  s=1; b=1
  resubmit signal || s=0
  resubmit bkgmc  || b=0
  sleep 60
  q=$(inqueue_cnt)
  log "pass $PASS: signal_resub=$s bkg_resub=$b in-queue-now=$q"
  if [ "$s" -eq 0 ] && [ "$b" -eq 0 ] && [ "$q" -eq 0 ]; then
    log "no more missing jobs -> all complete."
    break
  fi
  wait_drain
done

# ---- collect ----
log "running collector ..."
set +u; source "$CONDA_SH"; conda activate hza_ana; set -u
export PYTHONPATH="$REPO"; cd "$REPO"
python scripts/5_collect_cutflow.py >>"$OLOG" 2>&1
log "collector done. eveto keys check (signal 2022preEE):"
python - >>"$OLOG" 2>&1 <<'PY'
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
python scripts/plot_evetosigniVmA.py >>"$OLOG" 2>&1
log "plotter exit=$?"
ls -la "$PLOTREPO/plots/evetosigniVmA/" >>"$OLOG" 2>&1
log "orchestrator DONE. log=$OLOG"
