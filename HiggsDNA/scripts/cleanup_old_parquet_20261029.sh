#!/usr/bin/env bash
# cleanup_old_parquet_20261029.sh — 刪除使用者 2026-10-29 指定的 HZa 舊 parquet 產物（/eos/project）
# 理由：photon ID SF 換成自產 Z->ee 量測、muon SF 誤差改 RMS/sqrt(N) 後要重產 Sig/Bkg MC 與 control sample；
#   使用者勾選可刪：sfunc_smoke、sfunc、sfunc_y2024、fsrfix(Sig_MC)、zee_zmmg_ctrl(6 月 control)、parquet_eveto。
#   保留：parquet_DNA_tmp_sfunc_ord（目前 signal，重產後排列順序的參照）、fsrfix_fpo1（目前 bkg/data）等其他目錄。
# eos rm -r 對物件數大的目錄會拒絕但回 rc=0 → 每個目錄刪完都用 eos stat 驗證已消失；
#   沒消失就往下一層逐一刪，再刪外殼（最多遞迴 4 層）。
# 用法：bash cleanup_old_parquet_20261029.sh [--go]
set -uo pipefail
export EOS_MGM_URL=root://eosproject.cern.ch
P=/eos/project/h/htozg-dy-privatemc/pelai/HZa
TARGETS=(parquet_DNA_tmp_sfunc_smoke parquet_DNA_tmp_sfunc parquet_DNA_tmp_sfunc_y2024 parquet_DNA_tmp_fsrfix parquet_DNA_zee_zmmg_ctrl parquet_eveto)
GO=0; [ "${1:-}" = --go ] && GO=1

exists() { eos stat "$1" >/dev/null 2>&1; }
rmtree() {  # $1 path, $2 depth
  local d=$1 depth=$2
  eos rm -r "$d" > /dev/null 2>&1
  exists "$d" || return 0
  [ "$depth" -ge 4 ] && { echo "!! still exists at depth limit: $d"; return 1; }
  local c rc=0
  while read -r c; do [ -n "$c" ] && { rmtree "$d/$c" $((depth+1)) || rc=1; }; done < <(eos ls "$d" 2>/dev/null)
  eos rm -r "$d" > /dev/null 2>&1
  exists "$d" && { echo "!! failed: $d"; return 1; }
  return $rc
}

for t in "${TARGETS[@]}"; do
  exists "$P/$t" || { echo "already gone: $t"; continue; }
  echo "$t: $(eos du -s "$P/$t" 2>/dev/null | awk '{printf "%.1f GB", $1/1e9}')"
done
[ $GO = 1 ] || { echo "DRY-RUN; add --go to delete"; exit 0; }
bad=0
for t in "${TARGETS[@]}"; do
  exists "$P/$t" || continue
  echo "$(date +%T) deleting $t"
  rmtree "$P/$t" 0 || bad=$((bad+1))
  exists "$P/$t" && { echo "!! $t still exists"; bad=$((bad+1)); } || echo "$(date +%T) gone: $t"
done
echo "remaining under $P:"; eos ls "$P"
echo "failures $bad"; exit $bad
