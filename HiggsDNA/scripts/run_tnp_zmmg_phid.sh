#!/usr/bin/env bash
# =============================================================================
# run_tnp_zmmg_phid.sh — Z→μμγ photon ID TnP 生產（probe 不先要求 ID；ID 判定存成 probe_photon_pass_phid_custom_tight）
# =============================================================================
# 為什麼：原本 photon ID SF 用 EGM 的 Z→ee ntuple，量的 ID 和分析不同（chIso 絕對/相對、neutral iso/hcalPFClusterIso、
# H/E 的 EA 抄錯），且 10–20 GeV failing probe 背景佔 88–97%。改用 Z→μμγ：tagger 的 ID 與分析的
# photon_selections.select_resolved_photons 逐顆一致（studies/phid_zmmg_20261025/check_id_equivalence.py：2024 MC、2023 data 各 0 不一致）。
# config：metadata/za_{data,mc}_tnp_zmmg_phid_run3.json（probe_require_id=False、m_μμγ 窗 60–120 GeV）。
# CSEV 用的 za_*_tnp_zmmg_run3.json 與 tagger 預設行為不變。
#
# 用法（預設 dry-run，加 GO=1 才執行）：
#   MODE=data|mc  OUTDIR=...  YEARS=2024  SAMPLE_LIST=...  SHORT=1  BATCH_SYSTEM=local|condor  N_CORES=2  FPO=2  GO=1
# 預設輸出：/eos/cms/store/group/phys_susy/pelai/HZa/parquet_tnp_zmmg_phid/<MODE>
# log：HiggsDNA 寫在 OUTDIR 底下（job_* 目錄與 analysis_manager.pkl）；本腳本的 stdout 由呼叫端導向。
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
cd "${repo_dir}"

mode="${MODE:?MODE=data or mc}"
case "${mode}" in
  data) config="metadata/za_data_tnp_zmmg_phid_run3.json"; def_samples="Data_tnp_zmmg"; def_fpo=2 ;;
  mc)   config="metadata/za_mc_tnp_zmmg_phid_run3.json";   def_samples="DYJetsToLL,DYJetsTo2Mu,DYJetsToLL_MLM"; def_fpo=5 ;;
  *) echo "MODE must be data or mc"; exit 2 ;;
esac
outdir="${OUTDIR:-/eos/cms/store/group/phys_susy/pelai/HZa/parquet_tnp_zmmg_phid/${mode}}"
years="${YEARS:-2022preEE,2022postEE,2023preBPix,2023postBPix,2024}"
sample_list="${SAMPLE_LIST:-${def_samples}}"
fpo="${FPO:-${def_fpo}}"
n_cores="${N_CORES:-10}"
batch_system="${BATCH_SYSTEM:-condor}"
short="${SHORT:-0}"
clean_analysis_state="${CLEAN_ANALYSIS_STATE:-1}"
unretire_jobs="${UNRETIRE_JOBS:-0}"
merge_outputs="${MERGE_OUTPUTS:-1}"
# chunk size 1 會讓 poller hang（ref_hza_zee_zmmg_control），用 500
export HIGGSDNA_CONDOR_SUBMIT_CHUNK_SIZE="${CONDOR_SUBMIT_CHUNK_SIZE:-500}"
export HIGGSDNA_CONDOR_Q_TIMEOUT="${CONDOR_Q_TIMEOUT:-60}"

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate /eos/home-p/pelai/App/Conda/.conda/envs/hza_ana
set -u
export PYTHONNOUSERSITE=1
unset CPATH CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export PYTHONPATH="${repo_dir}"
export X509_USER_PROXY="${X509_USER_PROXY:-/afs/cern.ch/user/p/pelai/x509up_u$(id -u)}"

cmd=(python -s scripts/run_analysis.py --config "${config}" --log-level INFO --n_cores "${n_cores}"
     --output_dir "${outdir}" --batch_system "${batch_system}" --sample_list "${sample_list}" --years "${years}" --fpo "${fpo}"
     --reconfigure_jobs)
[[ "${unretire_jobs}" == "1" ]] && cmd+=(--unretire_jobs)
[[ "${merge_outputs}" == "1" ]] && cmd+=(--merge_outputs)
[[ "${short}" == "1" ]] && cmd+=(--short)

echo "[$(date '+%F %T')] mode=${mode} outdir=${outdir} years=${years} samples=${sample_list} fpo=${fpo} batch=${batch_system} short=${short}"
echo "cmd: ${cmd[*]}"
if [[ "${GO:-0}" != "1" ]]; then echo "DRY-RUN (set GO=1 to run)"; exit 0; fi
mkdir -p "${outdir}"
if [[ "${clean_analysis_state}" == "1" ]]; then
  rm -f "${outdir}/analysis_manager.pkl" "${outdir}/analysis_manager_temp.pkl"
fi
"${cmd[@]}"
