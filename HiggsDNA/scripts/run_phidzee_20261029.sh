#!/usr/bin/env bash
# =============================================================================
# run_phidzee_20261029.sh — HZa preselection 重產：Sig_MC + Bkg_MC + Z->ee/Z->mumugamma control sample
# =============================================================================
# 為什麼：photon ID SF 換成自產 Z->ee ntuple 的量測（分析實際用的 ID，使用者 10-29 目視全部擬合後定案），
#   muon SF 的擬合誤差改成 RMS/sqrt(N)。SF 與其 up/down 在 HiggsDNA 生產時寫進權重，所以 MC 要重產；
#   data 不受影響（使用者選 Sig + Bkg MC）。control sample（data + MC）一起重產，用來驗證新的 photon ID SF。
#   job 在執行當下從 AFS 讀 SF JSON：2024 的 photon ID SF 要等 DY_NLO 定案版重算完，所以用 YEARS 分兩批送。
#
# 用法：bash run_phidzee_20261029.sh <unit> <smoke|full> [--go]   （預設 dry-run）
#   YEARS=... 覆寫年份（預設五個 era）
#   unit：sig | bkg | bkg24x | ctrlmc | ctrlel | ctrlmu
#     sig     14 個 mA                                → parquet_DNA_tmp_phidzee/Sig_MC
#     bkg     DYGto2LG_10to100、DYJetsToLL、DYJetsTo2E/2Mu/2Tau（同 fsrfix_fpo1）→ parquet_DNA_tmp_phidzee/Bkg_MC
#     bkg24x  2024 的 TTG、TTGG_Run3、TTto2L2Nu、ZGG（同 fsrfix_extraBkg） → parquet_DNA_tmp_phidzee/Bkg_MC_extra2024
#     ctrlmc  control MC（DY；腳本預設清單的 DYGto2LG_10to50 在 2022preEE 沒有檔案→0 job 當掉，改用 bkg 那組）                          → parquet_DNA_zee_zmmg_ctrl/mc
#     ctrlel  control data，EGamma PD（Z->ee 區）        → parquet_DNA_zee_zmmg_ctrl/data_egamma
#     ctrlmu  control data，Muon PD（Z->mumugamma 區）   → parquet_DNA_zee_zmmg_ctrl/data_muon
#   smoke = SHORT=1（每個 sample x era 一個 job），輸出到 <base>_smoke。
# 設定（doc/HZa/hza_preselection_production_pitfalls.md 第 18 節）：FPO=1、HZA_STAGE_INPUTS=1、HZA_BMODE_PACK、
#   記憶體一律 3000（1 核 slot 池大得多；超過 3 GB 會被湊成多核、幾乎排不到。少數 OOM 的 job 個別 qedit 調高後 release，qedit 不觸發湊核）、CONDOR_SUBMIT_CHUNK_SIZE=500、systemd-run --user。
# log：HiggsDNA/signal_rerun_logs/phidzee_20261029/<unit>_<mode>_driver.log；condor log 在 HiggsDNA/eos_logs/。
set -o pipefail
UNITNAME="${1:-}"; MODE="${2:-}"; GO=0; [ "${3:-}" = "--go" ] && GO=1
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA
P=/eos/project/h/htozg-dy-privatemc/pelai/HZa
MA=mA_M1,mA_M2,mA_M3,mA_M4,mA_M5,mA_M6,mA_M7,mA_M8,mA_M9,mA_M10,mA_M15,mA_M20,mA_M25,mA_M30
ALLY=2022preEE,2022postEE,2023preBPix,2023postBPix,2024
case "$UNITNAME" in
  sig)    SCRIPT=run_ana_signal.sh;                 SAMPLES=$MA; Y=${YEARS:-$ALLY}; BASE=parquet_DNA_tmp_phidzee; SUB=Sig_MC; MEM=3000 ;;
  bkg)    SCRIPT=run_ana_bkgmc.sh;                  SAMPLES=DYGto2LG_10to100,DYJetsToLL,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau; Y=${YEARS:-$ALLY}; BASE=parquet_DNA_tmp_phidzee; SUB=Bkg_MC; MEM=3000 ;;
  bkg24x) SCRIPT=run_ana_bkgmc.sh;                  SAMPLES=TTG,TTGG_Run3,TTto2L2Nu,ZGG; Y=${YEARS:-2024}; BASE=parquet_DNA_tmp_phidzee; SUB=Bkg_MC_extra2024; MEM=3000 ;;
  ctrlmc) SCRIPT=run_ana_control_zee_zmmg.sh;       SAMPLES=DYGto2LG_10to100,DYJetsToLL,DYJetsTo2E,DYJetsTo2Mu,DYJetsTo2Tau; Y=${YEARS:-$ALLY}; BASE=parquet_DNA_zee_zmmg_ctrl; SUB=mc; MEM=3000 ;;
  ctrlel) SCRIPT=run_ana_data_control_zee_zmmg.sh;  SAMPLES=Data; Y=${YEARS:-$ALLY}; BASE=parquet_DNA_zee_zmmg_ctrl; SUB=data_egamma; MEM=3000 ;;
  ctrlmu) SCRIPT=run_ana_data_control_zee_zmmg.sh;  SAMPLES=Data_tnp_zmmg; Y=${YEARS:-$ALLY}; BASE=parquet_DNA_zee_zmmg_ctrl; SUB=data_muon; MEM=3000 ;;
  *) echo "usage: $0 sig|bkg|bkg24x|ctrlmc|ctrlel|ctrlmu smoke|full [--go]"; exit 1 ;;
esac
case "$MODE" in smoke) SHORT=1; BASE=${BASE}_smoke ;; full) SHORT=0 ;; *) echo "mode must be smoke|full"; exit 1 ;; esac
YTAG=$(echo "$Y" | tr ',' '\n' | sed 's/20//' | tr -d '\n')
OUTDIR=$P/$BASE/$SUB
LOGDIR=$REPO/signal_rerun_logs/phidzee_20261029
UNIT=hza-phidzee-$UNITNAME-$MODE-$YTAG

PROXY_SRC=/tmp/x509up_u$(id -u)
left=$(voms-proxy-info -file "$PROXY_SRC" -timeleft 2>/dev/null || echo 0)
[ "${left:-0}" -gt 86400 ] || { echo "!! proxy 剩 ${left}s（< 1 天），請先 voms-proxy-init"; exit 2; }

SAMPLE_ENV=""; [ -n "$SAMPLES" ] && SAMPLE_ENV="SAMPLE_LIST=$SAMPLES"
if [ "$GO" -eq 0 ]; then
  echo "DRY-RUN unit=$UNITNAME mode=$MODE script=$SCRIPT years=$Y ${SAMPLE_ENV}"
  echo "  OUTDIR=$OUTDIR  FPO=1 MEM=$MEM SHORT=$SHORT  systemd unit=$UNIT"
  echo "  log: $LOGDIR/${UNITNAME}_${MODE}_${YTAG}_driver.log  (proxy left ${left}s)"
  exit 0
fi
systemctl --user is-active "$UNIT" >/dev/null 2>&1 && { echo "!! unit $UNIT already active"; exit 3; }
mkdir -p "$LOGDIR" "$OUTDIR"
cp "$PROXY_SRC" "$REPO/x509up_u$(id -u)" && chmod 600 "$REPO/x509up_u$(id -u)"

RUNNER="$LOGDIR/${UNITNAME}_${MODE}_${YTAG}_runner.sh"
cat > "$RUNNER" <<EOF
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
cd $REPO
export PYTHONPATH=$REPO
export X509_USER_PROXY=$REPO/x509up_u$(id -u)
export HZA_BMODE_PACK=root://eoscms.cern.ch//eos/cms/store/group/phys_susy/pelai/App/hza_ana_pack.tar.gz
export HZA_STAGE_INPUTS=1
export HIGGSDNA_CONDOR_REQ_MEMORY=$MEM
echo "[\$(date '+%F %T')] start unit=$UNITNAME mode=$MODE years=$Y python=\$(which python) CONDA_PREFIX=\$CONDA_PREFIX"
OUTDIR=$OUTDIR $SAMPLE_ENV YEARS=$Y FPO=1 CONDOR_REQ_MEMORY=$MEM SHORT=$SHORT CLEAN_ANALYSIS_STATE=1 \\
  UNRETIRE_JOBS=1 MERGE_OUTPUTS=1 RECONFIGURE_JOBS=0 CONDOR_SUBMIT_CHUNK_SIZE=500 N_CORES=6 bash scripts/$SCRIPT
echo "[\$(date '+%F %T')] driver exit=\$?"
EOF
chmod +x "$RUNNER"
systemd-run --user --unit="$UNIT" --collect /bin/bash -c "exec /bin/bash $RUNNER >> $LOGDIR/${UNITNAME}_${MODE}_${YTAG}_driver.log 2>&1"
sleep 2
systemctl --user is-active "$UNIT"
