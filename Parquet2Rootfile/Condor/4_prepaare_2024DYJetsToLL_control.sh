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
DYJetsToLL_path=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control/DYJetsToLL
DYJetsTo2E_path=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control/DYJetsTo2E
DYJetsTo2Mu_path=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control/DYJetsTo2Mu
DYJetsTo2Tau_path=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_control/DYJetsTo2Tau

DYJetsToLL_path="$(strip_cr "$DYJetsToLL_path")"
DYJetsTo2E_path="$(strip_cr "$DYJetsTo2E_path")"
DYJetsTo2Mu_path="$(strip_cr "$DYJetsTo2Mu_path")"
DYJetsTo2Tau_path="$(strip_cr "$DYJetsTo2Tau_path")"

# 2024：把三個 channel hadd 成 DYJetsToLL/2024.root（不再使用 /<mA>/ 子目錄）
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