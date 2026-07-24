#!/usr/bin/env bash
# 等 HZa photon e-veto 研究的 signal+bkg condor jobs 跑完 → collect cutflow → plot。
# 只監看 eos_logs/Sig_MC 與 eos_logs/Bkg_MC 的 pelai jobs（不動使用者其他 job）。
set -uo pipefail

REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
PLOTREPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
LOGDIR="$REPO/eveto_submit_logs"
STAMP=$(date +%Y%m%d_%H%M%S)
WLOG="$LOGDIR/watch_${STAMP}.log"
mkdir -p "$LOGDIR"

log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$WLOG"; }

CONDA_SH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] || CONDA_SH="$(dirname "$(dirname "${CONDA_EXE:-$(command -v conda)}")")/etc/profile.d/conda.sh"

# active(idle+running) count of my eveto jobs
active_cnt(){ condor_q pelai -af JobStatus Cmd 2>/dev/null \
  | grep -E "eos_logs/(Sig_MC|Bkg_MC)/" | awk '$1==1||$1==2' | wc -l; }
held_cnt(){ condor_q pelai -af JobStatus Cmd 2>/dev/null \
  | grep -E "eos_logs/(Sig_MC|Bkg_MC)/" | awk '$1==5' | wc -l; }

log "watcher start; waiting for bkg submission to appear in queue..."
# 1) 等兩批都進 condor（至少 bkg 也出現，或已有大量 job）
for _ in $(seq 1 60); do
  nb=$(condor_q pelai -af Cmd 2>/dev/null | grep -cE "eos_logs/Bkg_MC/")
  ns=$(condor_q pelai -af Cmd 2>/dev/null | grep -cE "eos_logs/Sig_MC/")
  log "in-queue Sig=$ns Bkg=$nb"
  [ "$nb" -gt 0 ] && break
  sleep 60
done

# 2) 等 active jobs drain
STABLE=0
while :; do
  a=$(active_cnt); h=$(held_cnt)
  log "active(idle+running)=$a held=$h"
  if [ "$a" -eq 0 ]; then
    STABLE=$((STABLE+1))
    # 連續三次都 0 才確認（避免 submit 空窗 / condor 重排誤判）
    [ "$STABLE" -ge 3 ] && break
  else
    STABLE=0
  fi
  sleep 300
done
log "all active eveto jobs drained. held remaining=$(held_cnt)"

# 3) collect cutflow（hza_ana，純 pandas/numpy）
log "running collector (5_collect_cutflow.py) in hza_ana ..."
set +u; source "$CONDA_SH"; conda activate hza_ana; set -u
export PYTHONPATH="$REPO"
cd "$REPO"
python scripts/5_collect_cutflow.py >>"$WLOG" 2>&1
log "collector done. checking eveto keys in a signal JSON ..."
python - <<'PY' 2>>"$WLOG" | tee -a "$WLOG"
import json,glob
fs=sorted(glob.glob("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list/cutflow_Sig_MC_mA_M10_2018.json"))
fs=fs or sorted(glob.glob("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/cutflow_list/cutflow_Sig_MC_*_2022preEE.json"))
if fs:
    d=json.load(open(fs[0])); cf=d.get("cutflows",{})
    for k in ("zgammas_ph_eveto_w","zgammas_ph_noeveto_w"):
        c=cf.get(k,{}); print("  ",k,"allcuts=",c.get("all cuts") or c.get("event"))
PY
set +u; conda deactivate; set -u

# 4) plot（higgs-zg-plot_python，有 PyROOT）
log "running plotter (plot_evetosigniVmA.py) in higgs-zg-plot_python ..."
set +u; source "$CONDA_SH"; conda activate higgs-zg-plot_python; set -u
cd "$PLOTREPO"
python scripts/plot_evetosigniVmA.py >>"$WLOG" 2>&1
log "plotter exit=$?"
ls -la "$PLOTREPO/plots/evetosigniVmA/" 2>>"$WLOG" | tee -a "$WLOG"
log "watcher DONE. log: $WLOG"
