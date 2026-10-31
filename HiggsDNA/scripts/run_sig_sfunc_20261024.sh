#!/usr/bin/env bash
# =============================================================================
# run_sig_sfunc_20261024.sh — Run 3 signal 重產（photon ID 與 2024 electron SF 誤差修正後）
# =============================================================================
# 為什麼：SF 的 up/down 權重是 HiggsDNA 生產時寫進 parquet 的。TnP exportJson 原本只把
# nominal 擬合的統計誤差寫進 JSON，修正後（photon ID 全 era、electron 2024）要重產 signal
# 才會反映。同時帶入 09-27 之後更新過的 CSEV／electron／trigger SF。
# 選擇條件沒有改，所以新舊產物的事件必須逐筆相同，只有權重不同（由 compare 腳本驗證）。
#
# 用法：
#   bash run_sig_sfunc_20261024.sh smoke        # dry-run：印出小批測試（mA_M10 × 2023preBPix）的指令
#   bash run_sig_sfunc_20261024.sh smoke --go   # 實際送出小批測試
#   bash run_sig_sfunc_20261024.sh full  --go   # 全量 14 mA × 5 era（約 1104 jobs，fpo=1，memory 9000）
#   bash run_sig_sfunc_20261024.sh y2024 --go   # 只重產 2024（electron systdata 不含 statMC 之後）→ parquet_DNA_tmp_sfunc_y2024
#
# 輸出：smoke → /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_sfunc_smoke/Sig_MC
#       full  → /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_sfunc/Sig_MC
#       （舊產物 parquet_DNA_tmp_fsrfix/Sig_MC 不動，留作比對與退路）
# log ：HiggsDNA/signal_rerun_logs/sfunc_20261024/<mode>_driver.log；condor log 在 HiggsDNA/eos_logs/Sig_MC/
# 守護：driver 用 systemd-run --user 跑（lxplus 登出會清 process）。
set -o pipefail
MODE="${1:-}"; GO=0; [ "${2:-}" = "--go" ] && GO=1
case "$MODE" in
  smoke) SAMPLES=mA_M10; YEARS=2023preBPix; BASE=parquet_DNA_tmp_sfunc_smoke ;;
  full)  SAMPLES=mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30
         YEARS=2022preEE,2022postEE,2023preBPix,2023postBPix,2024; BASE=parquet_DNA_tmp_sfunc ;;
  # y2024：electron systdata 改成不含 statMC（10-24 22:5x）後，只有 2024 signal 受影響，單獨重產到新目錄
  y2024) SAMPLES=mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30
         YEARS=2024; BASE=parquet_DNA_tmp_sfunc_y2024 ;;
  *) echo "usage: $0 smoke|full|y2024 [--go]"; exit 1 ;;
esac

REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
OUTDIR=/eos/project/h/htozg-dy-privatemc/pelai/HZa/$BASE/Sig_MC
LOGDIR=$REPO/signal_rerun_logs/sfunc_20261024
UNIT=hza-sig-sfunc-$MODE

# proxy：repo 那份會被 HiggsDNA 帶進 job，要用有效的那份覆蓋（驗有效期，不是驗存在）
PROXY_SRC=/tmp/x509up_u$(id -u)
left=$(voms-proxy-info -file "$PROXY_SRC" -timeleft 2>/dev/null || echo 0)
[ "${left:-0}" -gt 86400 ] || { echo "!! proxy 剩 ${left}s（< 1 天），請先 voms-proxy-init"; exit 2; }

if [ "$GO" -eq 0 ]; then
  echo "DRY-RUN mode=$MODE samples=$SAMPLES years=$YEARS"
  echo "  OUTDIR=$OUTDIR"
  echo "  cp $PROXY_SRC $REPO/x509up_u$(id -u)"
  echo "  systemd-run --user --unit=$UNIT ... scripts/run_ana_signal.sh  (FPO=1 CONDOR_REQ_MEMORY=9000 CLEAN_ANALYSIS_STATE=1 HZA_BMODE_PACK=...)"
  echo "  log: $LOGDIR/${MODE}_driver.log"
  exit 0
fi

mkdir -p "$LOGDIR" "$OUTDIR"
cp "$PROXY_SRC" "$REPO/x509up_u$(id -u)" && chmod 600 "$REPO/x509up_u$(id -u)"

RUNNER="$LOGDIR/${MODE}_runner.sh"
cat > "$RUNNER" <<EOF
#!/usr/bin/env bash
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
cd $REPO
export PYTHONPATH=$REPO
export X509_USER_PROXY=$REPO/x509up_u$(id -u)
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
export HZA_STAGE_INPUTS=1
echo "[\$(date '+%F %T')] start mode=$MODE python=\$(which python) CONDA_PREFIX=\$CONDA_PREFIX"
OUTDIR=$OUTDIR SAMPLE_LIST=$SAMPLES YEARS=$YEARS FPO=1 CONDOR_REQ_MEMORY=9000 CLEAN_ANALYSIS_STATE=1 \\
  UNRETIRE_JOBS=1 MERGE_OUTPUTS=1 N_CORES=6 bash scripts/run_ana_signal.sh
echo "[\$(date '+%F %T')] driver exit=\$?"
EOF
chmod +x "$RUNNER"
systemd-run --user --unit="$UNIT" --collect /bin/bash -c "exec /bin/bash $RUNNER >> $LOGDIR/${MODE}_driver.log 2>&1"
sleep 2
systemctl --user is-active "$UNIT"
