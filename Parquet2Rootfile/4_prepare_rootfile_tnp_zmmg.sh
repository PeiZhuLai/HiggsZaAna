#!/usr/bin/env bash
set -Eeuo pipefail

# Prepare tnp_zmmg ROOT files produced by 3_run_P2Root_tnp_zmmg.sh.
# Input : /eos/home-p/pelai/HZa/root_P2Root/run3_tnp_zmmg_tmp
# Output: /eos/home-p/pelai/HZa/root_P2Root/run3_tnp_zmmg

strip_cr() { printf '%s' "$1" | tr -d $'\r'; }

INPUT_BASE="$(strip_cr /eos/home-p/pelai/HZa/root_P2Root/run3_tnp_zmmg_tmp)"
OUTPUT_BASE="$(strip_cr /eos/home-p/pelai/HZa/root_P2Root/run3_tnp_zmmg)"
YEARS=( 2022preEE 2022postEE 2023preBPix 2023postBPix 2024 )
DYJETS_LEGACY_YEARS=( 2022preEE 2022postEE 2023preBPix 2023postBPix )

AVAILABLE_MODULES=( data dyjets dyjets-mlm tt-dyjets )

print_usage() {
  cat <<'EOF_USAGE'
Usage:
  4_prepare_rootfile_tnp_zmmg.sh [--only <module1,module2,...>] [--list]

Modules:
  data        - Copy Data_tnp_zmmg yearly files and build run3.root
  dyjets      - Build DYJetsToLL yearly files; 2024 comes from DYJetsTo2Mu/2024.root
  dyjets-mlm  - Copy DYJetsToLL_MLM yearly files and build run3.root
  tt-dyjets   - hadd TT yearly files with prepared DYJetsToLL yearly files

Examples:
  4_prepare_rootfile_tnp_zmmg.sh
  4_prepare_rootfile_tnp_zmmg.sh --only dyjets,tt-dyjets
  4_prepare_rootfile_tnp_zmmg.sh --list
EOF_USAGE
}

selected_modules=()
if [[ $# -gt 0 ]]; then
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --only|-m)
        [[ $# -ge 2 ]] || { echo "ERROR: --only needs a module list"; exit 2; }
        IFS=',' read -r -a selected_modules <<< "${2// /,}"
        shift 2
        ;;
      --list)
        printf '%s\n' "${AVAILABLE_MODULES[@]}"
        exit 0
        ;;
      -h|--help)
        print_usage
        exit 0
        ;;
      *)
        echo "Unknown argument: $1"
        print_usage
        exit 2
        ;;
    esac
  done
fi

if [[ ${#selected_modules[@]} -eq 0 ]]; then
  selected_modules=( "${AVAILABLE_MODULES[@]}" )
fi

require_file() {
  local file
  file="$(strip_cr "$1")"
  if [[ ! -s "$file" ]]; then
    echo "ERROR: missing or empty file: $file"
    ls -l "$(dirname "$file")" || true
    exit 1
  fi
}

prepare_dir() {
  local dir
  dir="$(strip_cr "$1")"
  mkdir -p "$dir"
}

copy_root() {
  local source target
  source="$(strip_cr "$1")"
  target="$(strip_cr "$2")"
  require_file "$source"
  prepare_dir "$(dirname "$target")"
  echo "cp -f $source $target"
  cp -f "$source" "$target"
}

hadd_root() {
  local output
  output="$(strip_cr "$1")"
  shift

  local inputs=()
  local input
  for input in "$@"; do
    input="$(strip_cr "$input")"
    require_file "$input"
    inputs+=( "$input" )
  done

  prepare_dir "$(dirname "$output")"
  rm -f "$output"
  echo "hadd -f $output ${inputs[*]}"
  pushd "$(dirname "$output")" >/dev/null
  hadd -f "$(basename "$output")" "${inputs[@]}"
  popd >/dev/null
}

hadd_run3() {
  local sample_dir
  sample_dir="$(strip_cr "$1")"

  local inputs=()
  local year
  for year in "${YEARS[@]}"; do
    year="$(strip_cr "$year")"
    inputs+=( "$sample_dir/${year}.root" )
  done

  hadd_root "$sample_dir/run3.root" "${inputs[@]}"
}

prepare_data() {
  local out_dir="$OUTPUT_BASE/Data_tnp_zmmg"
  prepare_dir "$out_dir"

  local year
  for year in "${YEARS[@]}"; do
    year="$(strip_cr "$year")"
    copy_root "$INPUT_BASE/Data_tnp_zmmg/${year}.root" "$out_dir/${year}.root"
  done

  hadd_run3 "$out_dir"
}

prepare_dyjets() {
  local out_dir="$OUTPUT_BASE/DYJetsToLL"
  prepare_dir "$out_dir"

  local year
  for year in "${DYJETS_LEGACY_YEARS[@]}"; do
    year="$(strip_cr "$year")"
    copy_root "$INPUT_BASE/DYJetsToLL/${year}.root" "$out_dir/${year}.root"
  done

  copy_root "$INPUT_BASE/DYJetsTo2Mu/2024.root" "$out_dir/2024.root"
  hadd_run3 "$out_dir"
}

prepare_dyjets_mlm() {
  local out_dir="$OUTPUT_BASE/DYJetsToLL_MLM"
  prepare_dir "$out_dir"

  local year
  for year in "${YEARS[@]}"; do
    year="$(strip_cr "$year")"
    copy_root "$INPUT_BASE/DYJetsToLL_MLM/${year}.root" "$out_dir/${year}.root"
  done

  hadd_run3 "$out_dir"
}

prepare_tt_dyjets() {
  local dyjets_dir="$OUTPUT_BASE/DYJetsToLL"
  local out_dir="$OUTPUT_BASE/TT_DYJetsToLL"
  prepare_dir "$out_dir"

  local year
  for year in "${YEARS[@]}"; do
    year="$(strip_cr "$year")"
    hadd_root "$out_dir/${year}.root" \
      "$INPUT_BASE/TT/${year}.root" \
      "$dyjets_dir/${year}.root"
  done

  hadd_run3 "$out_dir"
}

for module in "${selected_modules[@]}"; do
  case "$module" in
    data)       prepare_data ;;
    dyjets)     prepare_dyjets ;;
    dyjets-mlm) prepare_dyjets_mlm ;;
    tt-dyjets)  prepare_tt_dyjets ;;
    *)
      echo "ERROR: unknown module '$module'"
      print_usage
      exit 2
      ;;
  esac
done
