#!/usr/bin/env bash
set -Eeuo pipefail

# 移除可能混入的 \r（CRLF）字元，避免 ROOT 將路徑截斷為 "/xxx.root"
strip_cr() { printf '%s' "$1" | tr -d $'\r'; }

# 新增：模組清單與使用說明、參數解析
# AVAILABLE_MODULES=( dygto2lg dyjets all-bkg data sig )
AVAILABLE_MODULES=( dyjets )

print_usage() {
  cat <<'EOF'
Usage:
  2_prepare_rootfile.sh [--only <module1,module2,...>] [--list]

Modules:
  dygto2lg   - Prepare annoying DYGto2LG
  dyjets     - Prepare DYJetsToLL
  all-bkg    - Add run3.root all background (需要 dygto2lg 與 dyjets 已完成)
  data       - Prepare Data
  sig        - Prepare Sig

Examples:
  4_prepare_rootfile.sh                 # 跑全部模組
  4_prepare_rootfile.sh --only dyjets   # 只跑 DYJetsToLL
  4_prepare_rootfile.sh -m dygto2lg,dyjets,all-bkg
  4_prepare_rootfile.sh --list
EOF
}

selected_modules=()
if [[ $# -gt 0 ]]; then
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --only|-m)
        [[ $# -ge 2 ]] || { echo "ERROR: --only 需要參數"; exit 2; }
        IFS=',' read -r -a selected_modules <<< "${2// /,}"
        shift 2
        ;;
      --list)
        printf '%s\n' "${AVAILABLE_MODULES[@]}"
        exit 0
        ;;
      -h|--help)
        print_usage; exit 0 ;;
      *)
        echo "Unknown argument: $1"
        print_usage
        exit 2
        ;;
    esac
  done
fi
# 預設跑全部模組
if [[ ${#selected_modules[@]} -eq 0 ]]; then
  selected_modules=( "${AVAILABLE_MODULES[@]}" )
fi

prepare_dyjets() {
#########################################################################
# # Prepare DYJetsToLL
#########################################################################
# 改用你提供的 run3_BDT 路徑
DYJetsToLL_path=/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/DYJetsToLL
DYJetsTo2E_path=/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/DYJetsTo2E
DYJetsTo2Mu_path=/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/DYJetsTo2Mu
DYJetsTo2Tau_path=/eos/home-p/pelai/HZa/root_P2Root/run3_BDT/DYJetsTo2Tau

DYJetsToLL_path="$(strip_cr "$DYJetsToLL_path")"
DYJetsTo2E_path="$(strip_cr "$DYJetsTo2E_path")"
DYJetsTo2Mu_path="$(strip_cr "$DYJetsTo2Mu_path")"
DYJetsTo2Tau_path="$(strip_cr "$DYJetsTo2Tau_path")"

# 新增：mA 1..30 各自把 2E/2Mu/2Tau 的 2024.root hadd 成 DYJetsToLL/<mA>/2024.root
mA_list=( \
  mA_M1 mA_M2 mA_M3 mA_M4 mA_M5 mA_M6 mA_M7 mA_M8 mA_M9 mA_M10 \
  mA_M11 mA_M12 mA_M13 mA_M14 mA_M15 mA_M16 mA_M17 mA_M18 mA_M19 mA_M20 \
  mA_M21 mA_M22 mA_M23 mA_M24 mA_M25 mA_M26 mA_M27 mA_M28 mA_M29 mA_M30 \
)

for ma in "${mA_list[@]}"; do
  ma="$(strip_cr "$ma")"

  in2e="$DYJetsTo2E_path/$ma/2024.root"
  in2mu="$DYJetsTo2Mu_path/$ma/2024.root"
  in2tau="$DYJetsTo2Tau_path/$ma/2024.root"
  outdir="$DYJetsToLL_path/$ma"
  out="$outdir/2024.root"

  # 檢查來源檔案存在且非空
  for f in "$in2e" "$in2mu" "$in2tau"; do
    if [ ! -s "$f" ]; then
      echo "ERROR: missing or empty file: $f"
      ls -l "$(dirname "$f")" || true
      exit 1
    fi
  done

  mkdir -p "$outdir"
  rm -f "$out"

  pushd "$outdir" >/dev/null
  echo "hadd $out $in2e $in2mu $in2tau"
  # 在 outdir 內輸出相對路徑檔名也可，但這裡直接用絕對路徑最不容易出錯
  hadd "2024.root" "$in2e" "$in2mu" "$in2tau"
  popd >/dev/null
done

# 2024：先把三個 channel hadd 成 DYJetsToLL/2024.root
for f in "$DYJetsTo2E_path/2024.root" "$DYJetsTo2Mu_path/2024.root" "$DYJetsTo2Tau_path/2024.root"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing or empty file: $f"
        ls -l "$(dirname "$f")" || true
        exit 1
    fi
done

rm -f "$DYJetsToLL_path/2024.root"
pushd "$DYJetsToLL_path" >/dev/null
echo "hadd $DYJetsToLL_path/2024.root $DYJetsTo2E_path/2024.root $DYJetsTo2Mu_path/2024.root $DYJetsTo2Tau_path/2024.root"
hadd "2024.root" "$DYJetsTo2E_path/2024.root" "$DYJetsTo2Mu_path/2024.root" "$DYJetsTo2Tau_path/2024.root"
popd >/dev/null
}


# 調度器：依指定模組順序執行
for m in "${selected_modules[@]}"; do
  case "$m" in
    dyjets)   prepare_dyjets ;;
    *)
      echo "ERROR: 未知模組 '$m'"; print_usage; exit 2 ;;
  esac
done