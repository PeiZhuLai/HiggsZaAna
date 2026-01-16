#!/usr/bin/env bash
set -Eeuo pipefail

# 移除可能混入的 \r（CRLF）字元，避免 ROOT 將路徑截斷為 "/xxx.root"
strip_cr() { printf '%s' "$1" | tr -d $'\r'; }

# 新增：模組清單與使用說明、參數解析
AVAILABLE_MODULES=( dygto2lg dyjets all-bkg data sig )
# AVAILABLE_MODULES=( sig )

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
  2_prepare_rootfile.sh                 # 跑全部模組
  2_prepare_rootfile.sh --only dyjets   # 只跑 DYJetsToLL
  2_prepare_rootfile.sh -m dygto2lg,dyjets,all-bkg
  2_prepare_rootfile.sh --list
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

# 將每個區塊包成函式（模組）

prepare_dygto2lg() {
#########################################################################
# # Prepare annoying DYGto2LG
#########################################################################
DYGto2LG_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG
DYGto2LG_10to50_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG_10to50
DYGto2LG_50to100_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG_50to100
DYGto2LG_10to100_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG_10to100

# 正規化路徑（移除 \r）
DYGto2LG_path="$(strip_cr "$DYGto2LG_path")"
DYGto2LG_10to50_path="$(strip_cr "$DYGto2LG_10to50_path")"
DYGto2LG_50to100_path="$(strip_cr "$DYGto2LG_50to100_path")"
DYGto2LG_10to100_path="$(strip_cr "$DYGto2LG_10to100_path")"

# hadd total.root file1.root file2.root
if [ -d "$DYGto2LG_path" ]; then
    echo "Directory exists: $DYGto2LG_path — removing it."
    rm -rf "$DYGto2LG_path"
else
    echo "Directory does not exist: $DYGto2LG_path — creating it."
fi
mkdir -p "$DYGto2LG_path"

echo "hadd $DYGto2LG_path/2022preEE.root $DYGto2LG_10to50_path/2022preEE.root $DYGto2LG_50to100_path/2022preEE.root"
hadd "$DYGto2LG_path/2022preEE.root" "$DYGto2LG_10to50_path/2022preEE.root" "$DYGto2LG_50to100_path/2022preEE.root"

echo "hadd $DYGto2LG_path/2022postEE.root $DYGto2LG_10to50_path/2022postEE.root $DYGto2LG_50to100_path/2022postEE.root"
hadd "$DYGto2LG_path/2022postEE.root" "$DYGto2LG_10to50_path/2022postEE.root" "$DYGto2LG_50to100_path/2022postEE.root"

# cp from.root to.root
echo "cp $DYGto2LG_10to100_path/2023preBPix.root $DYGto2LG_path/2023preBPix.root"
cp "$DYGto2LG_10to100_path/2023preBPix.root" "$DYGto2LG_path/2023preBPix.root"

echo "cp $DYGto2LG_10to100_path/2023postBPix.root $DYGto2LG_path/2023postBPix.root"
cp "$DYGto2LG_10to100_path/2023postBPix.root" "$DYGto2LG_path/2023postBPix.root"

echo "cp $DYGto2LG_10to100_path/2024.root $DYGto2LG_path/2024.root"
cp "$DYGto2LG_10to100_path/2024.root" "$DYGto2LG_path/2024.root"

# Sum four eras
years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

# 以陣列收集輸入檔，並對元素做 \r 清理
input_files=()
for year in "${years[@]}"; do
    y="$(strip_cr "$year")"
    input_files+=( "$DYGto2LG_path/${y}.root" )
done

# 檢查來源檔案存在且非空，方便快速定位問題
for f in "${input_files[@]}"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing or empty file: $f"
        ls -l "$(dirname "$f")" || true
        exit 1
    fi
done

# 在目標資料夾內用相對檔名執行 hadd，避免 EOS/xrootd 將絕對路徑誤解成 '/xxx.root'
base_inputs=()
for y in "${years[@]}"; do
    base_inputs+=( "$(strip_cr "$y").root" )
done

# 移除既有 run3.root 以避免舊檔干擾
rm -f "$DYGto2LG_path/run3.root"

pushd "$DYGto2LG_path" >/dev/null
echo "hadd $DYGto2LG_path/run3.root ${base_inputs[*]}"
hadd "run3.root" "${base_inputs[@]}"
popd >/dev/null
}

prepare_dyjets() {
#########################################################################
# # Prepare DYJetsToLL
#########################################################################
DYJetsToLL_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsToLL
DYJetsTo2E_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsTo2E
DYJetsTo2Mu_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsTo2Mu
DYJetsTo2Tau_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsTo2Tau

DYJetsToLL_path="$(strip_cr "$DYJetsToLL_path")"
DYJetsTo2E_path="$(strip_cr "$DYJetsTo2E_path")"
DYJetsTo2Mu_path="$(strip_cr "$DYJetsTo2Mu_path")"
DYJetsTo2Tau_path="$(strip_cr "$DYJetsTo2Tau_path")"

years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024 )

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

# 收集輸入檔，並對元素做 \r 清理
input_files=()
for year in "${years[@]}"; do
    y="$(strip_cr "$year")"
    input_files+=( "$DYJetsToLL_path/${y}.root" )
done

# 檢查來源檔案存在且非空
for f in "${input_files[@]}"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing or empty file: $f"
        ls -l "$(dirname "$f")" || true
        exit 1
    fi
done

# 在目標資料夾內用相對檔名執行 hadd
base_inputs=()
for y in "${years[@]}"; do
    base_inputs+=( "$(strip_cr "$y").root" )
done

# 移除既有 run3.root 以避免舊檔干擾（同 DYGto2LG）
rm -f "$DYJetsToLL_path/run3.root"

pushd "$DYJetsToLL_path" >/dev/null
echo "hadd $DYJetsToLL_path/run3.root ${base_inputs[*]}"
hadd "run3.root" "${base_inputs[@]}"
popd >/dev/null
}

prepare_all_bkg() {
#########################################################################
# Add run3.root all background
#########################################################################
# 依賴檢查：需要前兩個模組已產生 run3.root
if [ ! -s /eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsToLL/run3.root ] || \
   [ ! -s /eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG/run3.root ]; then
  echo "ERROR: 需要先完成模組 'dyjets' 與 'dygto2lg'，才可執行 'all-bkg'."
  exit 1
fi

DYJetsToLL_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYJetsToLL
DYGto2LG_path=/eos/home-p/pelai/HZa/root_P2Root/run3/DYGto2LG
Bkg_MC_path=/eos/home-p/pelai/HZa/root_P2Root/run3/All_Bkg
if [ -d "$Bkg_MC_path" ]; then
    echo "Directory exists: $Bkg_MC_path — removing it."
    rm -rf "$Bkg_MC_path"
else
    echo "Directory does not exist: $Bkg_MC_path — creating it."
fi
mkdir -p "$Bkg_MC_path"

echo "hadd $Bkg_MC_path/run3.root $DYJetsToLL_path/run3.root $DYGto2LG_path/run3.root"
hadd $Bkg_MC_path/run3.root $DYJetsToLL_path/run3.root $DYGto2LG_path/run3.root
}

prepare_data() {
#########################################################################
# # Prepare Data
#########################################################################
Data_path=/eos/home-p/pelai/HZa/root_P2Root/run3/Data
Data_path="$(strip_cr "$Data_path")"

years=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024)

# 收集輸入檔，並對元素做 \r 清理
input_files=()
for year in "${years[@]}"; do
    y="$(strip_cr "$year")"
    input_files+=( "$Data_path/${y}.root" )
done

# 檢查來源檔案存在且非空
for f in "${input_files[@]}"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing or empty file: $f"
        ls -l "$(dirname "$f")" || true
        exit 1
    fi
done

# 在目標資料夾內用相對檔名執行 hadd
base_inputs=()
for y in "${years[@]}"; do
    base_inputs+=( "$(strip_cr "$y").root" )
done

# 移除既有 run3.root 以避免舊檔干擾
rm -f "$Data_path/run3.root"

pushd "$Data_path" >/dev/null
echo "hadd $Data_path/run3.root ${base_inputs[*]}"
hadd "run3.root" "${base_inputs[@]}"
popd >/dev/null
}

prepare_sig() {
#########################################################################
# # Prepare Sig
#########################################################################
base_path=/eos/home-p/pelai/HZa/root_P2Root/run3
massList=( M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M15 M20 M25 M30 )
years=( 2024 )

# Add years into run3.root
for mass in "${massList[@]}"; do
    input_files=""
    for year in "${years[@]}"; do
        input_files+=" $base_path/mA_${mass}/${year}.root"
    done

    output_file="$base_path/mA_${mass}/run3.root"
    
    # Check if output file exists and remove it
    if [ -f "$output_file" ]; then
        echo "File exists: $output_file — removing it."
        rm -f "$output_file"
    fi
    
    echo "hadd $output_file $input_files"
    hadd $output_file $input_files
done

# Add run3.root in all ALP mass points
Sig_MC_path=/eos/home-p/pelai/HZa/root_P2Root/run3/All_Sig
if [ -d "$Sig_MC_path" ]; then
    echo "Directory exists: $Sig_MC_path — removing it."
    rm -rf "$Sig_MC_path"
else
    echo "Directory does not exist: $Sig_MC_path — creating it."
fi
mkdir -p "$Sig_MC_path"

input_files=""
for mass in "${massList[@]}"; do
    input_files+=" $base_path/mA_${mass}/run3.root"
done
output_file="$Sig_MC_path/run3.root"

# Check if output file exists and remove it
if [ -f "$output_file" ]; then
    echo "File exists: $output_file — removing it."
    rm -f "$output_file"
fi

echo "hadd $output_file $input_files"
hadd $output_file $input_files
}

# 調度器：依指定模組順序執行
for m in "${selected_modules[@]}"; do
  case "$m" in
    dygto2lg) prepare_dygto2lg ;;
    dyjets)   prepare_dyjets ;;
    all-bkg)  prepare_all_bkg ;;
    data)     prepare_data ;;
    sig)      prepare_sig ;;
    *)
      echo "ERROR: 未知模組 '$m'"; print_usage; exit 2 ;;
  esac
done