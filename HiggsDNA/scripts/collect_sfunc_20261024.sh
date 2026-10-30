#!/usr/bin/env bash
# =============================================================================
# collect_sfunc_20261024.sh — 只收 photon ID（全 era）與 electron 2024 重新輸出的 TnP JSON，再 merge
# =============================================================================
# 為什麼不跑 1_collect_custom_sf.sh：它也會從 TnP 輸出目錄重抓 muon、CSEV、electron 2025/26，
# 那些這次沒有重做；只收這次改過的 28 個檔，再跑 2_merge_custom_sf.py（輸入沒變的部分輸出必須
# 不變，事後用 git diff 驗證只有 unc_pass／systdata 類欄位改變）。
#
# 輸入：/eos/home-p/pelai/HZa/root_TnP/<flag>/<flag>.json
# 輸出：higgs_dna/systematics/data/<era>_UL/custom_SF_raw/<flag>.json → merge 成 hza_phid_*、hza_elid_2024、
#       hza_{dielleg12,dielleg23,sielleg30}trigger_2024、hza_eliso0p1/0p15_2024
# 預設 dry-run。實際執行加 --go。
set -euo pipefail
GO=0; [ "${1:-}" = "--go" ] && GO=1
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data
T=/eos/home-p/pelai/HZa/root_TnP
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/scripts

copy() {  # copy <era_dir> <flag>
  local src="$T/$2/$2.json" dst="$D/$1/custom_SF_raw/"
  [ -s "$src" ] || { echo "!! missing $src"; exit 3; }
  if [ "$GO" -eq 1 ]; then cp -p "$src" "$dst"; else echo "DRY-RUN cp $src $dst"; fi
}
for era in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do
  copy ${era}_UL hza_resolve_phid_${era}_sf
  copy ${era}_UL hza_resolve_phid_lowpt_${era}_sf
done
copy 2023postBPix_UL hza_resolve_phid_2023postBPixHole_sf
copy 2023postBPix_UL hza_resolve_phid_lowpt_2023postBPixHole_sf
for m in elid_nongap elid_gap elid_nongap_highpT elid_nongap_lowpT \
         elminiIso0p1_nongap elminiIso0p1_gap elminiIso0p1_nongap_highpT \
         elminiIso0p15_nongap elminiIso0p15_gap elminiIso0p15_nongap_highpT \
         dielleg12trigger_gap dielleg12trigger_nongap dielleg23trigger_gap dielleg23trigger_nongap \
         sielleg30trigger_gap sielleg30trigger_nongap; do
  copy 2024_UL hza_${m}_2024_sf
done
if [ "$GO" -eq 1 ]; then
  python3 "$S/2_merge_custom_sf.py" --electron-eras 2024 2025 2026
else
  echo "DRY-RUN python3 $S/2_merge_custom_sf.py --electron-eras 2024 2025 2026"
fi
